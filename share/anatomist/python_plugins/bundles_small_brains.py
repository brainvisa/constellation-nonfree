# -*- coding: utf-8 -*-

from __future__ import print_function

# Anatomist module
from __future__ import absolute_import
import anatomist.direct.api as anatomist

# BrainVisa module
from soma import aims

from soma.qt_gui.qt_backend import QtGui, QtCore

# Python system modules
import os
import sip
import math
import six

############################################################################
#                           Selection Action
############################################################################


class SmallBrainSelectionAction(anatomist.cpp.Action):
    def __init__(self):
        anatomist.cpp.Action.__init__(self)
        self.scaling = 0.1
        self.distance = 20  # mm
        self.display_all = True  # display all small brains on main view
        self._displayed_vertices = {}

    def name(self):
        return 'SmallBrainSelectionAction'

    def viewableAction(self):
        return False

    def addReducedObjectToView(self, aimsvertex, obj, nref, rot_center, tr, a,
                               window, diffuse_list):
        dup = obj.clone(True)  # shallow copy
        memo1 = self._displayed_vertices.setdefault(aimsvertex, {})
        memo = memo1.setdefault('objects', [])
        if dup is not None:
            a.registerObject(dup, False)
            dup.setReferential(nref)
            a.theProcessor().execute('SetMaterial', objects=[dup],
                                     diffuse=diffuse_list)
            window.registerObject(dup)
            memo.append([anatomist.cpp.rc_ptr_AObject(dup), rot_center, tr])
            a.releaseObject(dup)
            return dup
        else:
            if obj.objectTypeName(obj.type()) == 'TEXTURED SURF.':
                mesh = [x for x in obj
                        if x.objectTypeName(x.type()) == 'SURFACE'][0]
            else:
                mesh = obj
            mesh.setReferential(nref)
            a.theProcessor().execute('SetMaterial', objects=[mesh],
                                     diffuse=diffuse_list)
            mesh.notifyObservers()
            memo.append([anatomist.cpp.rc_ptr_AObject(obj), rot_center, tr])
            return obj

    def computeTransformation(self, graph, vertex):
        bb = (aims.Point3df(graph.boundingbox()[0][:3]),
              aims.Point3df(graph.boundingbox()[1][:3]))
        gcent = (bb[0] + bb[1])/2
        aimsgraph = graph.attributed()
        if 'normal' in vertex.attributed():
            normal = vertex.attributed()['normal']
            use_normal = True
        else:
            use_normal = False
        if not use_normal and 'lateralized_view' in aimsgraph:
            lateralized = True
            gcent_r = bb[0] + bb[1]
            gcent_r[0] *= 0.25
            gcent_r[1] /= 2
            gcent_r[2] /= 2
            gcent_l = bb[0] + bb[1]
            gcent_l[0] *= 0.75
            gcent_l[1] /= 2
            gcent_l[2] /= 2
        else:
            lateralized = False
        if 'center' in vertex.attributed():
            cent = aims.Point3df(vertex.attributed()['center'])
        else:
            bb = [aims.Point3df(x[:3]) for x in vertex.boundingbox()]
            if bb:
                cent = (bb[0] + bb[1])/2
            else:
                cent = aims.Point3df(128, 128, 64)
        if not use_normal:
            if lateralized:
                if (cent - gcent_l).norm2() > (cent - gcent_r).norm2():
                    gcent = gcent_r
                else:
                    gcent = gcent_l
            vdir = cent - gcent
            if vdir.norm2() != 0:
                vdir.normalize()
        else:
            vdir = normal
        m = aims.AffineTransformation3d()
        distance = self.distance
        if 'small_brains_distance' in aimsgraph:
            distance = aimsgraph['small_brains_distance']
        rot_center = cent + vdir * distance
        m.setTranslation(rot_center)
        scaling = self.scaling
        if 'small_brains_scaling' in aimsgraph:
            scaling = aimsgraph['small_brains_scaling']

        m *= aims.AffineTransformation3d([scaling, 0, 0, 0,
                                          0, scaling, 0, 0,
                                          0, 0, scaling, 0])
        m2 = aims.AffineTransformation3d()
        m2.setTranslation(-gcent)
        m *= m2
        return m, rot_center

    def computeRefAndTransformation(self, graph, vertex):
        m, rot_center = self.computeTransformation(graph, vertex)
        nref = anatomist.cpp.Referential()
        tr = anatomist.cpp.Transformation(nref, graph.getReferential())
        tr.motion().fromMatrix(m.toMatrix())
        sip.transferto(tr, None)  # don't delete tr
        tr.registerTrans()
        return nref, tr, rot_center

    def cleanup(self, keep=None):
        self.hide(keep=keep)
        a = anatomist.cpp.Anatomist()
        if hasattr(self, 'secondaryView'):
            secondview = self.secondaryView
            if secondview:
                objects = secondview.objects
                secondview.removeObjects(objects)
        refs = []
        for vertex, values in six.iteritems(self._displayed_vertices):
            if keep and vertex in keep:
                continue
            objs = values.get('objects', [])
            ref = values.get('referential', None)
            if ref is not None:
                refs.append(ref)
            for robj in objs:
                obj, rot_center, tr = robj
                graph = None
                parents = list(obj.parents())
                while parents:
                    p = parents.pop()
                    if p.type() == anatomist.cpp.AObject.GRAPH:
                        graph = p
                        break
                    parents += p.parents()
                if graph:
                    if obj.objectTypeName(obj.type()) == 'TEXTURED SURF.':
                        mesh = [x for x in obj.get()
                                if x.objectTypeName(x.type()) == 'SURFACE'][0]
                    else:
                        mesh = obj
                    mesh.setReferential(graph.getReferential())
        if not keep:
            self._displayed_vertices = {}
        else:
            self._displayed_vertices = dict(
                [(vertex, value)
                 for vertex, value in six.iteritems(self._displayed_vertices)
                 if vertex in keep])
        if refs:
            a.theProcessor().execute('DeleteElement', elements=refs)

    def hide(self, keep=None):
        a = anatomist.Anatomist()
        window = a.AWindow(a, self.view().aWindow())
        for vertex, values in six.iteritems(self._displayed_vertices):
            if keep and vertex in keep:
                continue
            objs = values.get('objects', [])
            for obj in objs:
                window.unregisterObject(obj[0].get())
            anavertex = vertex['ana_object']
            if anavertex:
                window.unregisterObject(anavertex)

    def initSmallBrains(self):
        window = self.view().aWindow()
        vertexlist = set()
        for obj in window.Objects():
            if obj.type() == anatomist.cpp.AObject.GRAPH:
                vertexlist = vertexlist.union(
                    [x.attributed() for x in obj
                     if isinstance(x.attributed(), aims.Vertex)])
        self.displayObjects(vertexlist, mode='set')

    def getClickedObject(self, x, y):
        window = self.view().aWindow()
        obj = window.objectAtCursorPosition(x, y)
        if obj is None:
            return None
        parents = [obj]
        while parents:
            parent = parents.pop(0)
            if parent.type() == anatomist.cpp.AObject.GRAPHOBJECT:
                return parent
            parents += list(parent.parents())
        for vertex, objs in six.iteritems(self._displayed_vertices):
            if obj in [o[0]._get() for o in objs['objects']]:
                return vertex['ana_object']
        return None

    def smallBrainClick(self, x, y, globX, globY):
        self.selectClickedSmallBrain(x, y, mode='set')

    def smallBrainClickAdd(self, x, y, globX, globY):
        self.selectClickedSmallBrain(x, y, mode='add')

    def smallBrainClickToggle(self, x, y, globX, globY):
        self.selectClickedSmallBrain(x, y, mode='toggle')

    def selectClickedSmallBrain(self, x, y, mode):
        sel_object = self.getClickedObject(x, y)
        window = self.view().aWindow()
        vertexlist = set()
        if sel_object:
            if sel_object.type() == anatomist.cpp.AObject.GRAPHOBJECT:
                go = sel_object.attributed()
                if isinstance(go, aims.Vertex):
                    vertexlist.add(go)
        if len(vertexlist) == 0:
            for obj in window.Objects():
                if obj.type() == anatomist.cpp.AObject.GRAPH:
                    vertexlist = vertexlist.union(
                        [xx.attributed() for xx in obj
                         if isinstance(xx.attributed(), aims.Vertex)])
        self.displayObjects(vertexlist, mode=mode)

    def displayObjects(self, vertexlist, mode):
        self.hide(keep=vertexlist)
        obj_to_display = []
        a = anatomist.Anatomist()
        window = self.view().aWindow()
        for aimsvertex in vertexlist:
            vmemo = self._displayed_vertices.setdefault(aimsvertex, {})
            vertex = aimsvertex['ana_object']
            # obj_to_display.append(vertex)
            objs = vmemo.get('objects', [])
            if len(objs) != 0:
                # vertex already done and available: just reuse it
                print('vertex in cache.')
                new_objects = [obj[0].get() for obj in objs]
                if not self.display_all:
                    # if display_all this will be done after the loop
                    for obj in new_objects:
                        window.registerObject(obj)
                obj_to_display += new_objects
                continue
            graph = list(vertex.parents())[0]
            nref, tr, rot_center = self.computeRefAndTransformation(graph,
                                                                    vertex)
            global_mesh = None
            if 'global_mesh' in graph.graph():
                global_mesh = graph.graph()['global_mesh']
            vmemo['referential'] = nref
            diffuse_list = [1., 0., 0., 1.]
            for obj in vertex:
                # print('obj:', type(obj), obj)
                obj_to_display.append(
                    self.addReducedObjectToView(aimsvertex,
                                                obj,
                                                nref,
                                                rot_center,
                                                tr,
                                                a,
                                                window,
                                                diffuse_list))
            if global_mesh is not None:
                diffuse_list = [0.8, 0.8, 0.8, 0.2]
                rmesh = self.addReducedObjectToView(aimsvertex, global_mesh,
                                                    nref,
                                                    rot_center,
                                                    tr,
                                                    a,
                                                    window,
                                                    diffuse_list)
                a.execute('SetMaterial', objects=[rmesh],
                          selectable_mode='always_selectable')
            for edge in aimsvertex.edges():
                if 'ana_object' in edge:
                    aedge = edge['ana_object']
                    diffuse_list = [0, 0, 1., 0.5]
                    for obj in aedge:
                        obj_to_display.append(self.addReducedObjectToView(
                            aimsvertex, obj, nref, rot_center, tr, a, window,
                            diffuse_list))

        if self.display_all:
            for vertex, values in six.iteritems(self._displayed_vertices):
                window.registerObject(vertex['ana_object'])
                for objd in values.get('objects', []):
                    window.registerObject(objd[0].get())

        if hasattr(self, 'secondaryView'):
            secondview = self.secondaryView
            if secondview:
                objects = secondview.objects
            if mode == 'set':
                secondview.removeObjects(objects)
                secondview.addObjects(obj_to_display)
            elif mode == 'add':
                secondview.addObjects(obj_to_display)
            elif mode == 'remove':
                secondview.removeObjects(obj_to_display)
            elif mode == 'toggle':
                to_add = []
                to_remove = []
                for obj in obj_to_display:
                    if secondview.hasObject(obj):
                        to_remove.append(obj)
                    else:
                        to_add.append(obj)
                if to_remove:
                    secondview.removeObjects(to_remove)
                if to_add:
                    secondview.addObjects(to_add)


class SmallBrainsRotationAction(anatomist.cpp.TrackOblique):
    def name(self):
        return 'SmallBrainsRotationAction'

    def beginTrackball(self, x, y, globalX, globalY):
        super(SmallBrainsRotationAction, self).beginTrackball(x,
                                                              y,
                                                              globalX,
                                                              globalY )
        self.tr_dict = {}
        action = self.view().controlSwitch().getAction(
            'SmallBrainSelectionAction')
        for vert_objs in six.itervalues(action._displayed_vertices):
            stored = vert_objs.get('objects', [])
            for obj, rot_center, tr in stored:
                if tr not in self.tr_dict:
                    self.tr_dict[tr] = aims.AffineTransformation3d(
                        tr.motion()), rot_center

    def moveTrackball(self, x, y, globalX, globalY):
        rot = self.rotation(x, y)
        v = rot.vector()
        v[3] *= -1
        rot.setVector(v)
        tr = None
        for tr, (begin_tr_AffineTransformation3d, rot_center) \
                in six.iteritems(self.tr_dict):
            rot_motion = aims.AffineTransformation3d()
            rot_motion.setTranslation(-rot_center)
            rot_center_translation = aims.AffineTransformation3d()
            rot_center_translation.setTranslation(rot_center)
            rot_motion \
                = (rot_center_translation *
                   aims.AffineTransformation3d(rot) *
                   rot_motion *
                   begin_tr_AffineTransformation3d)
            tr.motion().fromMatrix(rot_motion.toMatrix())
        if tr:
            tr.unregisterTrans()
            tr.registerTrans()
        self.view().aWindow().refreshNow()

    def cleanup(self):
        if hasattr(self, 'tr_dict'):
            del self.tr_dict


class SmallBrainsScaleAction(anatomist.cpp.Action):
    def name(self):
        return 'SmallBrainsScaleAction'

    def begin(self, x, y, globalX, globalY):
        self.tr_dict = {}
        self.y = y
        self.initscale = 1.
        action = self.view().controlSwitch().getAction(
            'SmallBrainSelectionAction')
        for vert_objs in six.itervalues(action._displayed_vertices):
            stored = vert_objs.get('objects', [])
            for obj, rot_center, tr in stored:
                if tr not in self.tr_dict:
                    self.tr_dict[tr] \
                        = aims.AffineTransformation3d(tr.motion()), rot_center

    def moveScale(self, x, y, globalX, globalY):
        ydiff = self.y - y
        zfac = math.exp(0.01 * ydiff)
        scale = self.initscale * zfac
        tr = None
        for tr, (begin_tr_AffineTransformation3d, rot_center) \
                in six.iteritems(self.tr_dict):
            scl_trans = aims.AffineTransformation3d()
            scl_trans.rotation().setValue(scale, 0, 0)
            scl_trans.rotation().setValue(scale, 1, 1)
            scl_trans.rotation().setValue(scale, 2, 2)
            rot_center_translation = aims.AffineTransformation3d()
            rot_center_translation.setTranslation(rot_center)
            scl_trans = rot_center_translation * scl_trans \
                * rot_center_translation.inverse() \
                * begin_tr_AffineTransformation3d
            tr.motion().fromMatrix(scl_trans.toMatrix())
        if tr:
            tr.unregisterTrans()
            tr.registerTrans()
        self.view().aWindow().refreshNow()

    def endScale(self, x, y, globalX, globalY):
        self.cleanup()

    def cleanup(self):
        if hasattr(self, 'tr_dict'):
            del self.tr_dict


class SmallBrainsTranslateAction(anatomist.cpp.Action):
    def name(self):
        return 'SmallBrainsTranslateAction'

    def begin(self, x, y, globalX, globalY):
        self.y = y
        action = self.view().controlSwitch().getAction(
            'SmallBrainSelectionAction')
        # what is the initial distance ? If specified in graphs, it should
        # be taken there, but what if several graphs have different distances
        # otherwise, default to the action distance
        self.initial_distance = action.distance
        graphs = [obj for obj in self.view().aWindow().Objects()
                  if obj.objectTypeName(obj.type()) == 'GRAPH']
        for graph in graphs:
            if 'small_brains_distance' in graph.attributed():
                self.initial_distance \
                    = graph.attributed()['small_brains_distance']
                break  # take just the first for now.

    def moveTranslate(self, x, y, globalX, globalY):
        # Warning: this action modifies /sets the graphs
        # 'small_brains_distance' attribute.
        trans_model = aims.Point3df()
        bbmin = self.view().windowBoundingMin()
        bbmax = self.view().windowBoundingMax()
        bby = bbmax[1] - bbmin[1]
        ydiff = float(self.y - y) / self.view().height() / self.view().zoom() \
            * bby
        action = self.view().controlSwitch().getAction(
            'SmallBrainSelectionAction')
        distance = self.initial_distance + ydiff
        tr = None
        done = {}
        for vertex, vert_objs in six.iteritems(action._displayed_vertices):
            ana_vertex = vertex['ana_object']
            graph = list(ana_vertex.parents())[0]
            graph.attributed()['small_brains_distance'] = distance
            stored = vert_objs.get('objects', [])
            for robj in stored:
                obj, rot_center, tr = robj
                new_rot_center = done.get(tr)
                if new_rot_center is None:
                    aims_tr, new_rot_center = action.computeTransformation(
                        graph, ana_vertex)
                    tr.motion().fromMatrix(aims_tr.toMatrix())
                    done[tr] = new_rot_center
                robj[1] = new_rot_center
        if tr:
            tr.unregisterTrans()
            tr.registerTrans()
        self.view().aWindow().refreshNow()

    def endTranslate(self, x, y, globalX, globalY):
        pass


class SmallBrainsControl(anatomist.cpp.Control3D):
    '''Control the view of a set of small glyph brains (textured brains).
    The control works with a graph representation of the regions. The graph
    should have a vertex for each region, each containing the representation of
    the region (mesh, textured mesh, bundles...)

    A reduced view of each region is displayed at the position of each region.
    An optional secondary view may be used to display a zoomed version of one
    or several small brains.
    Relations attached to a vertex are also displayed in the reduced view for
    this region.

    Interactions:
    - click on a small brain -> display it on the secondary view
    - click on the background -> display all reduced brains
    - shift + right mouse -> rotate all small views
    - ctrl + right mouse -> scale small brains
    - shift + ctrl + right mouse -> change the distance of the small brain from
      the center of the region

    Regions centers: either using the 'center' attribute of the graph vertices
    (if specified), or the their bounding boxes centers.

    Views are displayed at a distance of these centers:
    1. either moved along the normal, if vertices have the attribute 'normal'
    2. or moved along a line (graph_center - region_center)
       the center is the graph bounding box center.
    3. or moved along a line (hemisphere_center - region_center)
       the graph should have the attribute 'lateralized_view'.
       The left and right centers are determined by the graph bounding box:
       at 1/4 and 3/4 of the X axis, and middle of Y and Z.

    The distance is given either by the graph attribute 'small_brains_distance'
    or the SmallBrainSelectionAction action variable 'distance'.

    The reduction factor for small brains is given eitehr by the graph
    attribute 'small_brains_scaling', or the SmallBrainSelectionAction action
    variable 'scaling'.

    If the graph contains a global mesh in an attribute 'global_mesh', then
    this mesh is also displayed (reduced) in each small brain view (useful if
    the regions contain only fibers).
    '''
    def __init__(self):
        super(SmallBrainsControl, self).__init__(5792, 'SmallBrainsControl')

    def description(self):
        return '''<html><div align="center"><b>Small brains</b></div>
<p>This is a very specialized control, working on a graph built externally (BrainVisa viewers), containing a brain mesh, regions, and other things. It displays small copies of the brain mesh close to regions, each possibly having its own texture or fiber bundles attached. Selected small brain(s) can be displayed in a secondary window.
</p>
<p>
<table>
<tr><td>Left btn:</td><td>Select small brain, or display all if clicking on the background</td></tr>
<tr><td>+&&lt;shift&&gt;:</td><td>add small brain to 2nd view</td></tr>
<tr><td>+&&lt;ctrl&&gt;:</td><td>toggle small brain for 2nd view</td></tr>
<tr><td>&&lt;shift&&gt; + right btn:</td><td>rotate small brains</td></tr>
<tr><td>&&lt;ctrl&&gt; + right btn:</td><td>scale small brains</td></tr>
<tr><td>&&lt;shift&&gt; + &&lt;ctrl&&gt; + right btn:</td><td>change distance between region centers and small brains</td></tr>
</table>
</p>
</html>'''

    def eventAutoSubscription(self, pool):
        # plug parent control actions
        super(SmallBrainsControl, self).eventAutoSubscription(pool)
        # unplug the left mouse button action (normally used for linked cursor)
        # so that we can reuse the left button
        self.mouseLongEventUnsubscribe(
            QtCore.Qt.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier)
        self.mouseLongEventUnsubscribe(
            QtCore.Qt.RightButton, QtCore.Qt.ControlModifier)
        # now plug our new actions
        self.mousePressButtonEventSubscribe(
            QtCore.Qt.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier,
            pool.action('SmallBrainSelectionAction').smallBrainClick)
        self.mousePressButtonEventSubscribe(
            QtCore.Qt.LeftButton, QtCore.Qt.ShiftModifier,
            pool.action('SmallBrainSelectionAction').smallBrainClickAdd)
        self.mousePressButtonEventSubscribe(
            QtCore.Qt.LeftButton, QtCore.Qt.ControlModifier,
            pool.action('SmallBrainSelectionAction').smallBrainClickToggle)
        self.mouseLongEventSubscribe(
            QtCore.Qt.RightButton, QtCore.Qt.ShiftModifier,
            pool.action('SmallBrainsRotationAction').beginTrackball,
            pool.action('SmallBrainsRotationAction').moveTrackball,
            pool.action('SmallBrainsRotationAction').endTrackball, True)
        self.mouseLongEventSubscribe(
            QtCore.Qt.RightButton, QtCore.Qt.ControlModifier,
            pool.action('SmallBrainsScaleAction').begin,
            pool.action('SmallBrainsScaleAction').moveScale,
            pool.action('SmallBrainsScaleAction').endScale, True)
        self.mouseLongEventSubscribe(
            QtCore.Qt.RightButton,
            QtCore.Qt.ControlModifier | QtCore.Qt.ShiftModifier,
            pool.action('SmallBrainsTranslateAction').begin,
            pool.action('SmallBrainsTranslateAction').moveTranslate,
            pool.action('SmallBrainsTranslateAction').endTranslate, True)

    def doAlsoOnSelect(self, actionpool):
        super(SmallBrainsControl, self).doAlsoOnSelect(actionpool)
        actionpool.action('SmallBrainSelectionAction').initSmallBrains()

    def doAlsoOnDeselect(self, actionpool):
        super(SmallBrainsControl, self).doAlsoOnDeselect(actionpool)
        actionpool.action("SmallBrainsRotationAction").cleanup()
        actionpool.action("SmallBrainSelectionAction").cleanup()
        actionpool.action("SmallBrainsScaleAction").cleanup()

iconname = __file__
if iconname.endswith('.pyc') or iconname.endswith('.pyo'):
    iconname = iconname[:-1]
iconname = os.path.join(
    os.path.dirname(os.path.realpath(iconname)), 'bundles_small_brains.png')
pix = QtGui.QPixmap(iconname)
anatomist.cpp.IconDictionary.instance().addIcon('SmallBrainsControl', pix)
ad = anatomist.cpp.ActionDictionary.instance()
ad.addAction('SmallBrainSelectionAction', SmallBrainSelectionAction)
ad.addAction('SmallBrainsRotationAction', SmallBrainsRotationAction)
ad.addAction('SmallBrainsScaleAction', SmallBrainsScaleAction)
ad.addAction('SmallBrainsTranslateAction', SmallBrainsTranslateAction)
cd = anatomist.cpp.ControlDictionary.instance()
cd.addControl('SmallBrainsControl', SmallBrainsControl, 3792)
cm = anatomist.cpp.ControlManager.instance()
cm.addControl('QAGLWidget3D', '', 'SmallBrainsControl')
