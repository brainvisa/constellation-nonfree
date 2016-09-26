
from __future__ import print_function
import os, sys
import anatomist.direct.api as ana
from soma.qt_gui.qt_backend import QtGui, QtCore
from soma.qt_gui.qt_backend import init_matplotlib_backend
init_matplotlib_backend()
from matplotlib import pyplot
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

    def __init__(self, meshes, clusters, measurements, seed_gyri=[],
                 parent=None, flags=QtCore.Qt.WindowFlags(0)):
        '''Clusters inspector widget

        Parameters
        ----------
        meshes: list of aims.AimsTimeSurface objects
        clusters: list of aims.TimeTexture_S16 objects
        measurements: dict {int: pandas array}
        seed_gyri: list of aims.TimeTexture_S16 objects (optional)
        '''

        super(ClustersInspectorWidget, self).__init__(parent=parent,
                                                      flags=flags)
        if len(meshes) != len(clusters):
            raise ValueError('meshes and clusters numbers do not match')
        if len(seed_gyri) != 0 and len(seed_gyri) != len(meshes):
            raise ValueError('meshes and seed_gyri numbers do not match')

        a = ana.Anatomist('-b')
        self.meshes = [a.toAObject(mesh) for mesh in meshes]
        self.aims_clusters = clusters
        self.clusters = [a.toAObject(aims.TimeTexture('S16'))
                         for c in clusters]
        self.measurements = measurements
        self.seed_gyri = [a.toAObject(gyri) for gyri in seed_gyri]
        self.viewing_column = 0
        self.curves_columns = range(measurements[0].shape[1])

        # 3D views area
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

        # info dock
        info_dock = QtGui.QDockWidget()
        info_dock.setObjectName('info_dock')
        info_dock.setWindowTitle('Info')
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, info_dock)
        info = QtGui.QTextBrowser()
        info_dock.setWidget(info)
        self.info = info
        self.print_cluster_info(0, 1)

        # measurements table dock
        table_dock = QtGui.QDockWidget()
        table_dock.setObjectName('table_dock')
        table_dock.setWindowTitle('Table')
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

        # matplotlib curves dock
        curves_dock = QtGui.QDockWidget()
        curves_dock.setObjectName('curves_dock')
        curves_dock.setWindowTitle('Curves')
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, curves_dock)
        self.curves_fig = pyplot.figure()
        self.curves_widget = \
            pyplot._pylab_helpers.Gcf.get_fig_manager(
                self.curves_fig.number).window
        curves_dock.setWidget(self.curves_widget)
        toolbar = self.curves_widget.findChild(QtGui.QToolBar)
        toolbar.addSeparator()
        sel_cols = QtGui.QAction('Select plots', toolbar)
        toolbar.addAction(sel_cols)
        sel_cols.triggered.connect(self.select_curves_columns)

        # matrix view dock
        matrix_dock = QtGui.QDockWidget()
        matrix_dock.setObjectName('matrix_dock')
        matrix_dock.setWindowTitle('Matrix')
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, matrix_dock)
        matrix_dock.setWidget(QtGui.QLabel(
            'Here will be the matrix view.<br/>Soon.'))

        # fibers histogram dock
        fibers_histo_dock = QtGui.QDockWidget()
        fibers_histo_dock.setObjectName('fibers_histo_dock')
        fibers_histo_dock.setWindowTitle('Fibers length histogram')
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, fibers_histo_dock)
        self.fibers_histo_fig = pyplot.figure()
        self.fibers_histo_widget = \
            pyplot._pylab_helpers.Gcf.get_fig_manager(
                self.fibers_histo_fig.number).window
        fibers_histo_dock.setWidget(self.fibers_histo_widget)

        # cluster time evolution view
        cluster_time_dock = QtGui.QDockWidget()
        cluster_time_dock.setObjectName('cluster_time_dock')
        cluster_time_dock.setWindowTitle('Cluster evolution')
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, cluster_time_dock)
        self.cluster_time_fig = pyplot.figure()
        self.cluster_time_widget = \
            pyplot._pylab_helpers.Gcf.get_fig_manager(
                self.cluster_time_fig.number).window
        cluster_time_dock.setWidget(self.cluster_time_widget)

        # Anatomist views
        a = ana.Anatomist('-b')
        a.execute('GraphParams', show_tooltips=0)
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

        # setup a specific cluster time slider
        time_slider = clusters_win.findChild(QtGui.QSlider, 'sliderT')
        w = time_slider.parentWidget().parentWidget()
        slider_panel = QtGui.QWidget()
        lay = QtGui.QVBoxLayout()
        slider_panel.setLayout(lay)
        w.layout().addWidget(slider_panel)
        self.cluster_slider_label = QtGui.QLabel('2')
        self.cluster_slider = QtGui.QSlider()
        lay.addWidget(QtGui.QLabel('K:'))
        lay.addWidget(self.cluster_slider_label)
        lay.addWidget(self.cluster_slider)
        self.cluster_slider.setInvertedAppearance(True)
        self.cluster_slider.setInvertedControls(True)
        self.cluster_slider.setRange(2, len(self.measurements) + 1)
        self.cluster_slider.setPageStep(1)
        self.cluster_slider.valueChanged.connect(self.time_changed)

        self.clusters_fusions = []
        self.clusters_boundaries = []
        self.boundaries = []
        self.aims_measure_tex = []
        self.measure_tex = []
        self.measure_fusions = []

        self.update_clusters_texture()
        for mesh, clusters_tex, aims_clusters in zip(self.meshes,
                                                     self.clusters,
                                                     self.aims_clusters):
            clusters_tex.setPalette('Talairach')
            a.execute('TexturingParams', objects=[clusters_tex],
                      interpolation='rgb')
            clusters_fusion = a.fusionObjects(
                [mesh, clusters_tex], method='FusionTexSurfMethod')
            self.clusters_fusions.append(clusters_fusion)
            clusters_boundaries = a.toAObject(aims.AimsTimeSurface(2))
            self.clusters_boundaries.append(clusters_boundaries)
            clusters_boundaries.setMaterial(line_width=3.)

        clusters_win.addObjects(self.clusters_fusions)

        # remove referential button and views toolbar
        for win in (clusters_win, measures_win):
            win.findChild(QtGui.QPushButton).parent().hide()
            win.findChild(QtGui.QToolBar, 'mutations').hide()

        if len(self.seed_gyri) != 0:
            for mesh, seed_gyri_tex in zip(self.meshes, self.seed_gyri):
                boundaries \
                    = a.toAObject(aims.SurfaceManip.meshTextureBoundary(
                        mesh.surface(), seed_gyri_tex.toAimsObject(), -1))
                boundaries.setMaterial(line_width=3.)
                self.boundaries.append(boundaries)
            clusters_win.addObjects(self.boundaries)

        # build measurements texture and boundaries
        aims_tex = [aims.TimeTexture('FLOAT')] * len(self.meshes)
        self.aims_measure_tex = aims_tex
        self.measure_tex = [a.toAObject(tex) for tex in aims_tex]
        for measure_tex in self.measure_tex:
            measure_tex.setPalette('Yellow-red-fusion')
            measure_tex.setName(self.measurements[0].columns[self.viewing_column])
        self.make_measurements_texture()
        for mesh, measure_tex in zip(self.meshes, self.measure_tex):
            self.measure_fusions.append(a.fusionObjects(
                [mesh, measure_tex], method='FusionTexSurfMethod'))
        measures_win.addObjects(self.measure_fusions)
        # display colormap
        try:
            import paletteViewer
            # FIXME TODO have a common palette for all measure_tex
            paletteViewer.toggleShowPaletteForObject(self.measure_tex[0])
        except ImportError:
            pass
        self.update_clusters_boundaries()
        measures_win.addObjects(self.clusters_boundaries)

        install_controls()

        self.clusters_win.setControl('ClusterSelectionControl')

        # link time sliders on both views
        #time_slider.valueChanged.connect(self.time_changed)
        # Hide slider of 2nd window
        # hide both slider and label instead of hiding the parent panel
        # because the panel will be showed automatically by Window3D.Refresh()
        time_slider2 = measures_win.findChild(QtGui.QSlider, 'sliderT')
        time_slider2.hide()
        time_slider2.parentWidget().findChild(QtGui.QLabel).hide()

        # build and display table
        self.build_table(0)
        self.table.horizontalHeader().sectionDoubleClicked.connect(
            self.display_column)

        # build and display curves
        self.display_curves(0)


    def make_measurements_texture(self):
        col = self.viewing_column
        timestep = self.cluster_slider.value() - 2
        for atex, ctex, measure_tex in zip(self.aims_measure_tex,
                                           self.aims_clusters,
                                           self.measure_tex):
            ndata = len(ctex[0].data())
            arr = np.zeros((ndata,), dtype=np.float32)
            if len(ctex) > timestep:
                measurements = self.measurements[timestep]
                values = measurements[measurements.columns[col]]
                for label in range(len(values)):
                    value = values[label]
                    arr[ctex[timestep].data().arraydata() == label + 1] = value
            atex[0].data().assign(arr)
            measure_tex.setTexture(atex)
            measure_tex.setChanged()
            measure_tex.notifyObservers()


    def print_cluster_info(self, timestep, cluster):
        info_text = '''<h1>Cluster %d info:</h1>
Num of clusters (K): <b>%d</b><br/>
''' \
            % (cluster, timestep + 2)
        patch = self.patch_of_cluster(timestep, cluster)
        if patch is not None:
            if patch[1] is not None:
                info_text += 'In patch (gyrus): %d, label: %s<br/>' % patch
            else:
                info_text += 'In patch (gyrus): %d<br/>' % patch[0]
        self.info.setText(info_text)


    def patch_of_cluster(self, timestep, cluster):
        if len(self.seed_gyri) == 0:
            return None, None
        for aims_clusters, seed_gyri in zip(self.aims_clusters, self.seed_gyri):
            if len(aims_clusters) > timestep:
                vert = np.where(aims_clusters[timestep].arraydata() == cluster)
                if len(vert[0]) != 0:
                    patch_label = seed_gyri.toAimsObject()[0].data()[vert[0][0]]
                    print('patch_label:', patch_label)
                    hdr = seed_gyri.attributed()
                    if 'GIFTI_labels_table' in hdr:
                        label = hdr['GIFTI_labels_table'][patch_label]
                        print('gyrus label:', label)
                        return patch_label, label
                    return patch_label, None


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
        self.set_curve_cursor(timestep, label)
        self.select_table_row(timestep, label)
        self.update_cluster_time(timestep, label)
        self.clusters_win.statusBar().showMessage('Cluster: %d' % label)


    def update_clusters_boundaries(self):
        timestep = self.cluster_slider.value() - 2
        for mesh, aims_clusters, boundaries in zip(self.meshes,
                                                   self.aims_clusters,
                                                   self.clusters_boundaries):
            texture = aims.TimeTexture(dtype='S16')
            texture[0].data().assign(
                aims_clusters[min(timestep, len(aims_clusters) - 1)].data())
            clusters_boundaries \
                = aims.SurfaceManip.meshTextureBoundary(
                    mesh.surface(), texture, -1)
            boundaries.setSurface(clusters_boundaries)
            boundaries.notifyObservers()


    def update_clusters_texture(self):
        timestep = self.cluster_slider.value() - 2
        for clusters, aims_clusters in zip(self.clusters, self.aims_clusters):
            t = min(timestep, len(aims_clusters) - 1)
            texture = aims.TimeTexture('S16')
            texture[0].data().assign(aims_clusters[t].data())
            clusters.setTexture(texture)
            clusters.notifyObservers()


    def time_changed(self, value):
        self.cluster_slider_label.setText(str(value))
        timestep = value - 2
        self.build_table(timestep)
        self.display_curves(timestep)
        self.make_measurements_texture()
        self.update_clusters_texture()
        self.update_clusters_boundaries()


    def display_column(self, col):
        if self.viewing_column != col:
            name = self.measurements[0].columns[col]
            self.column_label.setText('displaying: <b>%s</b>' % name)
            self.viewing_column = col
            self.make_measurements_texture()
            try:
                import paletteViewer
                paletteViewer.toggleShowPaletteForObject(self.measure_tex[0])
                self.measure_tex[0].setName(name)
                #self.measure_tex[0].setChanged()
                #self.measure_tex[0].notifyObservers()
                paletteViewer.toggleShowPaletteForObject(self.measure_tex[0])
            except ImportError:
                pass

    def display_curves(self, timestep):
        if len(self.curves_fig.axes) == 0:
            # no axes yet. Create them
            bgcolor = self.palette().color(QtGui.QPalette.Base)
            axes = self.curves_fig.add_subplot(111, axisbg=str(bgcolor.name()))
        else:
            # use existing axes
            axes = self.curves_fig.axes[0]
            axes.clear()
        axes.set_xlabel('cluster')
        data = self.measurements[timestep]
        cols = data.columns[self.curves_columns]
        if len(cols) != 0:
            axes.set_xticks(range(data.shape[0] + 1))
            axes.set_xlim((0.8, data.shape[0] + 0.2))
            axes.plot(range(1, data.shape[0] + 1), data[cols], 'o-')
            axes.legend(list(cols))
        self.curves_fig.canvas.draw()


    def set_curve_cursor(self, timestep, label):
        if len(self.curves_fig.axes) != 0:
            axes = self.curves_fig.axes[0]
            if len(axes.lines) > len(self.curves_columns):
                # erase previous cursor
                axes.lines.remove(axes.lines[-1])
            lim = axes.get_ylim()
            lim = (lim[0] + lim[1] * 0.01, lim[1] * 0.99)
            axes.plot([label - 0.1, label - 0.1, label + 0.1, label + 0.1,
                       label - 0.1],
                      lim + lim[::-1] + (lim[0], ), 'r')
            self.curves_fig.canvas.draw()


    def select_curves_columns(self):
        dialog = QtGui.QDialog()
        dialog.setModal(True)
        lay = QtGui.QVBoxLayout()
        dialog.setLayout(lay)
        lwid = QtGui.QListWidget()
        lay.addWidget(lwid)
        lwid.addItems(self.measurements[0].columns)
        lwid.setSelectionMode(lwid.ExtendedSelection)
        for col in self.curves_columns:
            lwid.item(col).setSelected(True)
        ok = QtGui.QPushButton('OK')
        cancel = QtGui.QPushButton('Cancel')
        hlay = QtGui.QHBoxLayout()
        lay.addLayout(hlay)
        hlay.addStretch()
        hlay.addWidget(ok)
        hlay.addWidget(cancel)
        ok.pressed.connect(dialog.accept)
        cancel.pressed.connect(dialog.reject)
        res = dialog.exec_()
        if res:
            self.curves_columns = [item.row()
                                   for item in lwid.selectedIndexes()]
            timestep = int(round(self.clusters_win.getTime()))
            self.display_curves(timestep)


    def select_table_row(self, timestep, label):
        self.table.selectRow(label - 1)


    def update_cluster_time(self, timestep, label):
        if len(self.cluster_time_fig.axes) == 0:
            # no axes yet. Create them
            bgcolor = self.palette().color(QtGui.QPalette.Base)
            axes = self.cluster_time_fig.add_subplot(
                111, axisbg=str(bgcolor.name()))
        else:
            # use existing axes
            axes = self.cluster_time_fig.axes[0]
            axes.clear()
        axes.set_xlabel('cluster evolution')
        # TODO: put some data in the plot here...
        self.cluster_time_fig.canvas.draw()


def load_clusters_instpector_files(mesh_filenames, clusters_filenames,
                                   measurements_filenames,
                                   seed_gyri_filenames=[]):
    meshes = []
    clusters = []
    seed_gyri = []
    for mesh_filename in mesh_filenames:
        meshes.append(aims.read(mesh_filename))
    for clusters_filename in clusters_filenames:
        clusters.append(aims.read(clusters_filename))
    measurements = None
    for seed_gyri_filename in seed_gyri_filenames:
        seed_gyri.append(aims.read(seed_gyri_filename))

    if len(meshes) != len(clusters):
        raise ValueError('meshes and clusters numbers do not match')
    if len(seed_gyri) != 0 and len(seed_gyri) != len(meshes):
        raise ValueError('meshes and seed_gyri numbers do not match')
    return meshes, clusters, measurements, seed_gyri


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

    use_ex_num = 0

    if use_ex_num == 1:
        meshes, clusters, measurements, seed_gyri = load_clusters_instpector_files(
            #['/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.3.0/group_analysis/01to40/average_brain/averagebrain.white.mesh'],
            ['/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.1.0/001/surf/bh.r.aims.white.inflated.gii'],
            ['/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to40/connectivity_clustering/avg/fs01to40/lh.supramarginal/smooth3.0/avgSubject/01to40_avg_fs01to40_lh.supramarginal_avgSubject_clusteringTime.gii'],
            None,
            ['/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.3.0/group_analysis/01to40/average_brain/bh.annot.averagebrain.gii'])
    elif use_ex_num == 1:
        meshes, clusters, measurements, seed_gyri = load_clusters_instpector_files(
            ['/neurospin/population/HCP/S500-1/100307/T1w/fsaverage_LR32k/100307.L.inflated.32k_fs_LR.surf.gii', '/neurospin/population/HCP/S500-1/100307/T1w/fsaverage_LR32k/100307.R.inflated.32k_fs_LR.surf.gii'],
            ['/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to40/connectivity_clustering/avg/fs01to40/lh.supramarginal/smooth3.0/avgSubject/01to40_avg_fs01to40_lh.supramarginal_avgSubject_clusteringTime.gii', '/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to40/connectivity_clustering/avg/fs01to40/lh.supramarginal/smooth3.0/avgSubject/01to40_avg_fs01to40_lh.supramarginal_avgSubject_clusteringTime.gii'],
            None,
            ['/neurospin/population/HCP/S500-1/100307/MNINonLinear/fsaverage_LR32k/100307.L.aparc.32k_fs_LR.label.gii', '/neurospin/population/HCP/S500-1/100307/MNINonLinear/fsaverage_LR32k/100307.R.aparc.32k_fs_LR.label.gii'])
    # temp
    measurements = dict(
        (i, pandas.DataFrame(np.random.ranf((i + 2, 4)),
                             columns=('size', 'homogeneity', 'conn_density', 'other')))
        for i in range(max([len(clusters_tex) for clusters_tex in clusters])))
    cw = ClustersInspectorWidget(
        meshes, clusters, measurements=measurements, seed_gyri=seed_gyri)
    cw.show()

    if run_event_loop:
        QtGui.QApplication.instance().exec_()


