
from __future__ import print_function
import os, sys
import anatomist.direct.api as ana
from soma.qt_gui.qt_backend import QtGui, QtCore
from soma import aims
import numpy as np
import weakref
import pandas
import six

# need to initialize Anatomist before importing selection plugin
a = ana.Anatomist('-b')

import selection  # anatomist plugin

class ClusterSelectionAction(ana.cpp.Action):

    def __init__(self):
        super(ClusterSelectionAction, self).__init__()

    def name(self):
        return 'ClusterSelectionAction'

    def select_cluster(self, x, y, dx, dy):
        win = self.view().aWindow()
        obj = win.objectAtCursorPosition(x, y)
        if obj is None:
            return
        mesh = [m for m in obj if isinstance(m, ana.cpp.ASurface_3)][0]
        tex = [m for m in obj if isinstance(m, ana.cpp.ATexture)][0]
        poly = win.polygonAtCursorPosition(x, y, obj)
        polygon = mesh.surface().polygon()[poly]
        vertices = mesh.surface().vertex()
        pos = win.positionFromCursor(x, y)
        vert = [vertices[p] - pos for p in polygon]
        ivert = polygon[np.argmin([x.norm2() for x in vert])]
        a = ana.Anatomist('-b')
        timestep = win.getTime() / obj.TimeStep()
        aims_tex = a.AObject(a, tex).toAimsObject()
        timestep = max([x for x in aims_tex.keys() if x <= timestep])
        label = aims_tex[timestep][ivert]
        print('label:', label)
        if hasattr(win, 'cluster_window'):
            win.cluster_window.cluster_selected(timestep, label)


class ClusterSelectionControl(selection.SelectionControl):

    def __init__(self):
        super(ClusterSelectionControl, self).__init__(
            'ClusterSelectionControl')


    def eventAutoSubscription(self, pool):
        super(ClusterSelectionControl, self).eventAutoSubscription(pool)
        self.mousePressButtonEventUnsubscribe(
            QtCore.Qt.LeftButton, QtCore.Qt.NoModifier)
        self.mousePressButtonEventSubscribe(
            QtCore.Qt.LeftButton, QtCore.Qt.NoModifier,
            pool.action('ClusterSelectionAction' ).select_cluster)


class ClustersInspectorWidget(QtGui.QMainWindow):

    def __init__(self, mesh, clusters, measurements, seed_gyri=None,
                 parent=None, flags=QtCore.Qt.WindowFlags(0)):
        super(ClustersInspectorWidget, self).__init__(parent=parent,
                                                      flags=flags)
        self.mesh = mesh
        self.clusters = clusters
        self.aims_clusters = clusters.toAimsObject()
        self.measurements = measurements
        self.seed_gyri = seed_gyri
        self.viewing_column = 0

        main_w = QtGui.QSplitter(QtCore.Qt.Vertical)
        main_w.setObjectName('views')
        self.setCentralWidget(main_w)

        clusters_view = QtGui.QWidget()
        clusters_view.setObjectName('clusters_view')
        main_w.addWidget(clusters_view)
        clusters_view.setLayout(QtGui.QHBoxLayout())

        measures_view = QtGui.QWidget()
        measures_view.setObjectName('measures_view')
        main_w.addWidget(measures_view)
        measures_view.setLayout(QtGui.QHBoxLayout())

        info_dock = QtGui.QDockWidget()
        info_dock.setObjectName('info_dock')
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, info_dock)
        info = QtGui.QTextBrowser()
        info_dock.setWidget(info)
        self.info = info
        self.print_cluster_info(0, 1)

        table_dock = QtGui.QDockWidget()
        table_dock.setObjectName('table_dock')
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, table_dock)
        table_wid = QtGui.QWidget()
        table_dock.setWidget(table_wid)
        table_wid.setLayout(QtGui.QVBoxLayout())
        self.column_label = QtGui.QLabel(
            'displaying: <b>%s</b>'
            % measurements[0].columns[self.viewing_column])
        table_wid.layout().addWidget(self.column_label)
        table = QtGui.QTableWidget()
        table_wid.layout().addWidget(table)
        self.table = table

        a = ana.Anatomist('-b')
        clusters_win = a.createWindow('3D')
        clusters_view.layout().addWidget(clusters_win.getInternalRep())
        self.clusters_win = clusters_win
        clusters_win.getInternalRep().cluster_window = weakref.proxy(self)

        measures_win = a.createWindow('3D')
        measures_view.layout().addWidget(measures_win.getInternalRep())
        self.measures_win = measures_win

        a.execute('WindowConfig', windows=[clusters_win, measures_win],
                  linkedcursor_on_slider_change=1)
        a.execute('LinkWindows', windows=[clusters_win, measures_win],
                  group=153)

        clusters.setPalette('Talairach')
        a.execute('TexturingParams', objects=[clusters], interpolation='rgb')
        self.clusters_fusion = a.fusionObjects([mesh, clusters],
                                               method='FusionTexSurfMethod')
        self.clusters_boundaries \
            = a.toAObject(aims.SurfaceManip.meshTextureBoundary(
                mesh.surface(), self.aims_clusters, -1))
        self.clusters_boundaries.setMaterial(line_width=3.)
        measures_win.addObjects(self.clusters_boundaries)

        clusters_win.addObjects(self.clusters_fusion)

        if self.seed_gyri is not None:
            self.boundaries \
                = a.toAObject(aims.SurfaceManip.meshTextureBoundary(
                    mesh.surface(), self.seed_gyri.toAimsObject(), -1))
            self.boundaries.setMaterial(line_width=3.)
            clusters_win.addObjects(self.boundaries)

        aims_tex = aims.TimeTexture('FLOAT')
        self.aims_measure_tex = aims_tex
        self.measure_tex = a.toAObject(aims_tex)
        self.measure_tex.setPalette('Yellow-red-fusion')
        self.make_measurements_texture(0)
        self.measure_funsion = a.fusionObjects([mesh, self.measure_tex],
                                               method='FusionTexSurfMethod')
        measures_win.addObjects(self.measure_funsion)

        install_controls()

        self.clusters_win.setControl('ClusterSelectionControl')

        time_slider = clusters_win.findChild(QtGui.QSlider, 'sliderT')
        time_slider.valueChanged.connect(self.time_changed)

        self.build_table(0)
        self.table.horizontalHeader().sectionDoubleClicked.connect(
            self.display_column)


    def make_measurements_texture(self, col):
        atex = self.aims_measure_tex
        ctex = self.aims_clusters
        ndata = len(ctex[0].data())
        for i in ctex.keys():
            measurements = self.measurements[i]
            values = measurements[measurements.columns[col]]
            arr = np.zeros((ndata,), dtype=np.float32)
            for label in range(len(values)):
                value = values[label]
                arr[ctex[i].data().arraydata() == label + 1] = value
            atex[i].data().assign(arr)
        self.measure_tex.setTexture(atex)
        self.measure_tex.setChanged()
        self.measure_tex.notifyObservers()


    def print_cluster_info(self, timestep, cluster):
        info_text = '''<h1>Cluster %d info:</h1>
Nb of clusters: <b>%d</b><br/>
Timestep: <b>%d</b><br/>
''' \
            % (cluster, timestep + 2, timestep)
        self.info.setText(info_text)


    def build_table(self, timestep):
        self.table.clear()
        measurements = self.measurements[timestep]
        self.table.setColumnCount(measurements.shape[1])
        self.table.setRowCount(measurements.shape[0])
        self.table.setHorizontalHeaderLabels(measurements.columns)
        #hdr = self.table.horizontalHeader()
        for c in range(measurements.shape[1]):
            col = measurements[measurements.columns[c]]
            for r in range(measurements.shape[0]):
                self.table.setItem(r, c,
                                   QtGui.QTableWidgetItem(str(col[r])))


    def cluster_selected(self, timestep, label):
        self.print_cluster_info(timestep, label)


    def time_changed(self):
        timestep = self.clusters_win.getTime() / self.clusters.TimeStep()
        self.build_table(timestep)

    def display_column(self, col):
        if self.viewing_column != col:
            self.column_label.setText('displaying: <b>%s</b>'
                                      % self.measurements[0].columns[col])
            self.make_measurements_texture(col)
            self.viewing_column = col


def load_clusters_instpector_files(mesh_filename, clusters_filename,
                                   measurements_filenames,
                                   seed_gyri_filename=None):
    a = ana.Anatomist('-b')
    mesh = a.loadObject(mesh_filename)
    clusters = a.loadObject(clusters_filename)
    measurements = None
    seed_gyri = None
    if seed_gyri_filename is not None:
        seed_gyri = a.loadObject(seed_gyri_filename)

    return mesh, clusters, measurements, seed_gyri


def clear_controls():
    ad = ana.cpp.ActionDictionary.instance()
    ad.removeAction('ClusterSelectionAction')
    cd = ana.cpp.ControlDictionary.instance()
    cd.removeControl('ClusterSelectionControl')
    cm = ana.cpp.ControlManager.instance()
    cm.removeControl('QAGLWidget3D', '', 'ClusterSelectionControl')


def install_controls():
    clear_controls()

    #iconname = __file__
    #if iconname.endswith( '.pyc' ) or iconname.endswith( '.pyo' ):
        #iconname = iconname[:-1]
    #iconname = os.path.join(os.path.dirname(os.path.realpath(iconname)),
                            #'cluster_selection_icon.png')
    #pix = QtGui.QPixmap(iconname)
    pix = ana.cpp.IconDictionary.instance().getIconInstance('SelectionControl')
    ana.cpp.IconDictionary.instance().addIcon('ClusterSelectionControl',
                                                    pix )
    ad = ana.cpp.ActionDictionary.instance()
    ad.addAction('ClusterSelectionAction', ClusterSelectionAction)
    cd = ana.cpp.ControlDictionary.instance()
    cd.addControl('ClusterSelectionControl', ClusterSelectionControl, 100)
    cm = ana.cpp.ControlManager.instance()
    cm.addControl('QAGLWidget3D', '', 'ClusterSelectionControl')



if __name__ == '__main__':
    run_event_loop = False
    if QtGui.QApplication.instance() is None:
        qapp = QtGui.QApplication(['bloup'])
        run_event_loop = True

    mesh, clusters, measurements, seed_gyri = load_clusters_instpector_files(
        '/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.3.0/group_analysis/01to40/average_brain/averagebrain.white.mesh',
        '/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to40/connectivity_clustering/avg/fs01to40/lh.supramarginal/smooth3.0/avgSubject/01to40_avg_fs01to40_lh.supramarginal_avgSubject_clusteringTime.gii',
        None,
        '/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.3.0/group_analysis/01to40/average_brain/bh.annot.averagebrain.gii')
    # temp
    measurements = dict(
        (i, pandas.DataFrame(np.random.ranf((i + 2, 2)),
                             columns=('size', 'homogeneity')))
        for i in range(len(clusters.toAimsObject())))
    cw = ClustersInspectorWidget(
        mesh, clusters, measurements=measurements, seed_gyri=seed_gyri)
    cw.show()

    if run_event_loop:
        QtGui.QApplication.instance().exec_()


