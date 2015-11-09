# -*- coding: utf-8 -*-

# Anatomist module
import anatomist.direct.api as anatomist

# BrainVisa module
from soma import aims

from PyQt4 import QtGui, QtCore
import os, sip
import selection

class BundlesSelectionAction( selection.SelectionAction ):
  def __init__(self):
    selection.SelectionAction.__init__(self)
    self.scaling = 0.1
    self.decal = 20#mm

  def name( self ):
    return 'BundlesSelectionAction'

  def viewableAction( self ):
    return False

  def addReducedObjectToView(self, obj, nref, rot_center, tr, a, window, diffuse_list):
    dup = obj.clone( True ) # shallow copy
    print 'dup:', dup
    if dup is not None:
      a.registerObject( dup, False ) # False )
      dup.setReferential( nref )
      a.theProcessor().execute( 'SetMaterial', objects=[dup],
                diffuse=diffuse_list )
      window.registerObject( dup )
      self._stored.append( (anatomist.cpp.rc_ptr_AObject( dup ), rot_center, tr) )
      a.releaseObject( dup )
      return dup
    else:
      print 'no dup'
      if obj.objectTypeName( obj.type() ) == 'TEXTURED SURF.':
        mesh = [ x for x in obj \
          if x.objectTypeName( x.type() ) == 'SURFACE' ][0]
      else:
        mesh = obj
      print 'mesh:', mesh
      mesh.setReferential( nref )
      a.theProcessor().execute( 'SetMaterial', objects=[mesh],
                diffuse=diffuse_list )
      mesh.notifyObservers()
      self._stored.append( (anatomist.cpp.rc_ptr_AObject( obj ), rot_center, tr) )
      return obj

  def computeRefAndTransformation(self, graph, vertex):
    bb = graph.boundingbox()
    gcent = (bb[0] + bb[1])/2
    a = anatomist.cpp.Anatomist()
    nref = anatomist.cpp.Referential()
    if vertex.attributed().has_key( 'center' ):
      print 'using existing center'
      cent = aims.Point3df( vertex.attributed()[ 'center' ] )
      print 'center:', cent
    else:
      bb = vertex.boundingbox()
      cent = (bb[0] + bb[1])/2
    vdir = cent - gcent
    if vdir.norm2() != 0:
      vdir.normalize()
    tr = anatomist.cpp.Transformation( nref, graph.getReferential() )
    m = tr.motion()
    rot_center = cent + vdir * self.decal
    print 'rot_center:', rot_center
    m.setTranslation( rot_center )
    m *= aims.Motion( [ self.scaling, 0, 0, 0,
                        0, self.scaling, 0, 0,
                        0, 0, self.scaling, 0 ] )
    m2 = aims.Motion()
    m2.setTranslation( -gcent )
    m *= m2
    tr.motion().fromMatrix( m.toMatrix() )
    sip.transferto( tr, None ) # don't delete tr
    tr.registerTrans()
    return nref, tr, rot_center

  def releaseObjectsRefs( self ):
    if not hasattr( self, '_stored' ):
      self._stored = []
      return
    for obj, rot_center, tr in self._stored:
      graph = None
      parents = list( obj.parents() )
      print 'parents:', parents
      while parents:
        p = parents.pop()
        print 'test parent:', p
        if p.type() == anatomist.cpp.AObject.GRAPH:
          print 'parent graph found'
          graph = p
          break
        parents += p.parents()
      if graph:
        print 'reset ref of graph', graph
        if obj.objectTypeName( obj.type() ) == 'TEXTURED SURF.':
          print 'obj:', obj.name()
          print 'children:', [ x for x in obj.get() ]
          mesh = [ x for x in obj.get() \
            if x.objectTypeName( x.type() ) == 'SURFACE' ][0]
        else:
          mesh = obj
        print 'mesh:', mesh
        mesh.setReferential( graph.getReferential() )
      else: print 'graph not found'
    self._stored = []

  def edgeSelection( self ):
    try:
      recursion = getattr( self, '_recursion' )
      if recursion:
        # print "recursion"
        return
    except:
      pass
    selection.SelectionAction.edgeSelection( self )
    self._recursion = True
    selectedEdges_set = set()
    self.cleanup()
    sf = anatomist.cpp.SelectFactory.factory()
    sel = sf.selected()
    window = self.view().aWindow()
    group = window.Group()
    gsel = sel.get( group )
    vertexlist = set()
    edgeslist = set()
    a = anatomist.Anatomist()
    if gsel:
      for obj in gsel:
        print 'selected:', obj
        if obj.type() == anatomist.cpp.AObject.GRAPHOBJECT:
          print 'is a GRAPHOBJECT'
          go = obj.attributed()
          if isinstance( go, aims.Vertex ):
            print 'is a Vertex'
            vertexlist.add( go )
            #for edge in go.edges():
              #edgeslist.add(edge)
          else: print 'not a Vertex'
        else: print 'not a GRAPHOBJECT'
    else:
      print 'nothing selected'
      for obj in window.Objects():
        print 'obj:', obj
        if obj.type() == anatomist.cpp.AObject.GRAPH:
          print 'graph.'
          vertexlist = vertexlist.union( [ x.attributed() for x in obj \
            if isinstance( x.attributed(), aims.Vertex ) ] )
    print 'vertexlist:', len( vertexlist )
    print vertexlist
    self.releaseObjectsRefs()
    self._storedrefs = []
    obj_to_display = []
    for aimsvertex in vertexlist:
      vertex = aimsvertex['ana_object']
      obj_to_display.append( vertex )
      graph = list(vertex.parents())[0]
      nref, tr, rot_center = self.computeRefAndTransformation(graph, vertex)
      global_mesh = None
      if graph.graph().has_key('global_mesh'):
        global_mesh = graph.graph()['global_mesh']
      self._storedrefs.append( nref )
      diffuse_list = [ 1., 0., 0., 1. ]
      for obj in vertex:
        print 'obj:', type( obj ), obj
        obj_to_display.append(self.addReducedObjectToView(obj, nref, rot_center, tr, a, window, diffuse_list))
      if global_mesh is not None:
        diffuse_list = [ 0.8, 0.8, 0.8, 0.2 ]
        self.addReducedObjectToView(global_mesh, nref, rot_center, tr, a, window, diffuse_list)
      for edge in aimsvertex.edges():
      #for edge in edgeslist:
        if edge.has_key( 'ana_object' ):
          aedge = edge[ 'ana_object' ]
          diffuse_list = [ 0, 0, 1., 0.5 ]
          for obj in aedge:
            print 'add small edge object', obj.name()
            obj_to_display.append(self.addReducedObjectToView(obj, nref, rot_center, tr, a, window, diffuse_list))

    del self._recursion
    if hasattr( self, 'secondaryView' ):
      secondview = self.secondaryView
      print 'secondaryView:', secondview
      print 'objects:', [ type(x) for x in secondview.Objects() ]
      objects = secondview.objects
      print 'objects, 2:', objects
      secondview.removeObjects( objects )
      secondview.addObjects( obj_to_display )
    else:
      print 'no secondaryView in', self

  def cleanup( self ):
    self._recursion = True
    if hasattr( self, '_stored' ):
      self.releaseObjectsRefs()
    if hasattr( self, '_storedrefs' ):
      a = anatomist.cpp.Anatomist()
      a.theProcessor().execute( 'DeleteElement', elements=self._storedrefs )
      del self._storedrefs
    self._recursion = False


class BundlesRotationSelectionAction( anatomist.cpp.TrackOblique ):
  def name( self ):
    return 'BundlesRotationSelectionAction'

  def beginTrackball( self, x, y, globalX, globalY ):
    anatomist.cpp.TrackOblique.beginTrackball(self, x, y, globalX, globalY )
    self.tr_dict = {}
    action = self.view().controlSwitch().getAction( 'BundlesSelectionAction' )
    try:
      stored = action._stored
    except:
      return
    for obj, rot_center, tr in stored:
      if not self.tr_dict.has_key(tr):
        self.tr_dict[tr] = aims.Motion(tr.motion()), rot_center
  
  def moveTrackball( self, x, y, globalX, globalY ):
    rot = self.rotation(x, y)
    v = rot.vector()
    v[3] *= -1;
    rot.setVector( v )
    tr = None
    for tr, (begin_tr_Motion, rot_center) in self.tr_dict.iteritems():
      rot_motion = aims.Motion()
      rot_motion.setTranslation(-rot_center)
      rot_center_translation = aims.Motion()
      rot_center_translation.setTranslation(rot_center)
      rot_motion = rot_center_translation*aims.Motion(rot)*rot_motion*begin_tr_Motion
      tr.motion().fromMatrix( rot_motion.toMatrix() )
    if tr:
      tr.unregisterTrans()
      tr.registerTrans()
    self.view().aWindow().refreshNow()

  def cleanup( self ):
    if hasattr( self, 'tr_dict' ):
      del self.tr_dict


class BundlesSelectionControl( selection.SelectionControl ):
  def __init__( self ):
    selection.SelectionControl.__init__( self, 'BundlesSelectionControl' )

  def eventAutoSubscription( self, pool ):
    anatomist.cpp.Select3DControl.eventAutoSubscription( self, pool )
    self.selectionChangedEventSubscribe( pool.action( 'BundlesSelectionAction' ).edgeSelection )
    self.mouseLongEventUnsubscribe( QtCore.Qt.MidButton, QtCore.Qt.ShiftModifier )
    self.mouseLongEventSubscribe( \
    QtCore.Qt.MidButton, QtCore.Qt.ShiftModifier, pool.action( 'BundlesRotationSelectionAction' ).beginTrackball, pool.action( 'BundlesRotationSelectionAction' ).moveTrackball, pool.action( 'BundlesRotationSelectionAction' ).endTrackball , True )

  def doAlsoOnSelect( self, actionpool ):
    actionpool.action( 'BundlesSelectionAction' ).edgeSelection()

  def doAlsoOnDeselect( self, actionpool ):
    actionpool.action("BundlesRotationSelectionAction").cleanup()
    actionpool.action("BundlesSelectionAction").cleanup()
    selection.SelectionControl.doAlsoOnDeselect( self, actionpool )

iconname = __file__
if iconname.endswith( '.pyc' ) or iconname.endswith( '.pyo' ):
  iconname = iconname[:-1]
iconname = os.path.join(os.path.dirname(os.path.realpath(iconname)),'bundlesSelectionIcon.png')
pix = QtGui.QPixmap( iconname )
anatomist.cpp.IconDictionary.instance().addIcon( 'BundlesSelectionControl',
  pix )
ad = anatomist.cpp.ActionDictionary.instance()
ad.addAction( 'BundlesSelectionAction', BundlesSelectionAction )
ad.addAction( 'BundlesRotationSelectionAction', BundlesRotationSelectionAction )
cd = anatomist.cpp.ControlDictionary.instance()
cd.addControl( 'BundlesSelectionControl', BundlesSelectionControl, 30 )
cm = anatomist.cpp.ControlManager.instance()
cm.addControl( 'QAGLWidget3D', '', 'BundlesSelectionControl' )

