
from __future__ import print_function
from __future__ import absolute_import
import os, sys
import anatomist.direct.api as ana
from soma.qt_gui.qt_backend import QtGui, QtCore
from soma.qt_gui.qt_backend import init_matplotlib_backend
from six.moves import range
from six.moves import zip
init_matplotlib_backend()
from matplotlib import pyplot
from soma import aims
import numpy as np
import weakref
from soma.aims import pandas
import six
from constel.lib.utils import matrixtools


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
        if obj is None or (not isinstance(obj, ana.cpp.ASurface_3) 
                           and not isinstance(obj, ana.cpp.MObject)):
            return
        if isinstance(obj, ana.cpp.ASurface_3):
            mesh = obj
        else:
            meshes = [m for m in obj if isinstance(m, ana.cpp.ASurface_3)]
            if len(meshes) == 0:
                return
            mesh = meshes[0]
        tex = [m for m in obj if isinstance(m, ana.cpp.ATexture)][0]
        poly = win.polygonAtCursorPosition(x, y, obj)
        polygon = mesh.surface().polygon()[poly]
        vertices = mesh.surface().vertex()
        pos = win.positionFromCursor(x, y)
        vert = [vertices[p] - pos for p in polygon]
        ivert = polygon[np.argmin([x.norm2() for x in vert])]
        a = ana.Anatomist('-b')
        mesh = a.AObject(a, mesh)
        timestep = win.getTime() / obj.TimeStep()
        aims_tex = a.AObject(a, tex).toAimsObject()
        timestep = max([x for x in aims_tex.keys() if x <= timestep])
        label = aims_tex[timestep][ivert]
        if hasattr(win, 'cluster_window'):
            win.cluster_window.cluster_selected(timestep, label, mesh, ivert)


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
    '''Clusters inspector widget

    Parameters
    ----------
    meshes: list of aims.AimsTimeSurface objects
    clusters: list of aims.TimeTexture_S16 objects
    measurements: dict {int: pandas array}
    seed_gyri: list of aims.TimeTexture_S16 objects (optional)
    matrix: aims.SparseOrDenseMatrix (optional)
    parent: parent widget (optional)
    flags: Qt widget flags (optional)
    '''

    def __init__(self, meshes, clusters, measurements, seed_gyri=[],
                 matrix=None,
                 parent=None, flags=QtCore.Qt.WindowFlags(0)):

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
        self.aims_seed_gyri = seed_gyri
        self.seed_gyri = [a.toAObject(gyri) for gyri in seed_gyri]
        self.matrix = matrix
        self.clusters_matrix = {}
        self.viewing_column = 0
        self.curves_columns = list(range(list(measurements.values())[0].shape[1]))
        self.selected_cluster_boundaries = a.toAObject(aims.AimsTimeSurface(2))

        # 3D views area
        main_w = QtGui.QSplitter(QtCore.Qt.Vertical)
        main_w.setObjectName('views')
        self.setCentralWidget(main_w)

        self.create_info_dock()
        self.create_table_dock()
        self.create_curves_dock()
        self.create_anatomist_views()
        self.create_matrix_dock()
        self.create_fibers_histo_dock()
        self.create_clusters_evolution_dock()
        self.create_silhouette_dock()

        # build and display table
        self.build_table(list(measurements.keys())[0])
        self.table.horizontalHeader().sectionDoubleClicked.connect(
            self.display_column)

        # build and display curves
        self.display_curves(list(measurements.keys())[0])
        
        
    def __del__(self):
        print('##### DEL ClustersInspectorWidget #####')
        if self.curves_fig is not None:
            pyplot.close(self.curves_fig)
        if self.matrix_fig is not None:
            pyplot.close(self.matrix_fig)
        if self.fibers_histo_fig is not None:
            pyplot.close(self.fibers_histo_fig)
        if self.cluster_time_fig is not None:
            pyplot.close(self.cluster_time_fig)


    def create_info_dock(self):
        # info dock
        info_dock = QtGui.QDockWidget()
        info_dock.setObjectName('info_dock')
        info_dock.setWindowTitle('Info')
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, info_dock)
        info = QtGui.QTextBrowser()
        info_dock.setWidget(info)
        self.info = info
        self.print_cluster_info(0, 1)


    def create_table_dock(self):
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
            % list(self.measurements.values())[0].columns[self.viewing_column])
        table_wid.layout().addWidget(self.column_label)
        table = QtGui.QTableWidget()
        table_wid.layout().addWidget(table)
        self.table = table


    def create_curves_dock(self):
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


    def create_matrix_dock(self):
        # matrix view dock
        matrix_dock = QtGui.QDockWidget()
        matrix_dock.setObjectName('matrix_dock')
        matrix_dock.setWindowTitle('Matrix')
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, matrix_dock)
        self.matrix_fig = pyplot.figure()
        self.matrix_widget = \
            pyplot._pylab_helpers.Gcf.get_fig_manager(
                self.matrix_fig.number).window
        matrix_dock.setWidget(self.matrix_widget)
        self.display_matrix()
        toolbar = self.matrix_widget.findChild(QtGui.QToolBar)
        toolbar.hide()
        statusbar = self.matrix_widget.findChild(QtGui.QStatusBar)
        statusbar.hide()


    def create_fibers_histo_dock(self):
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


    def create_clusters_evolution_dock(self):
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


    def create_silhouette_dock(self):
        # silhouette dock
        silhouette_dock = QtGui.QDockWidget()
        silhouette_dock.setObjectName('silhouette_dock')
        silhouette_dock.setWindowTitle('Silhouette')
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, silhouette_dock)
        silhouette_dock.setWidget(QtGui.QLabel('Nice silhouette here. Cute.'))


    def create_anatomist_views(self):
        # Anatomist views
        a = ana.Anatomist('-b')
        main_w = self.centralWidget()

        a.execute('GraphParams', show_tooltips=0)

        clusters_view = QtGui.QWidget()
        clusters_view.setObjectName('clusters_view')
        main_w.addWidget(clusters_view)
        clusters_view.setLayout(QtGui.QHBoxLayout())

        clusters_win = a.createWindow('3D')
        clusters_view.layout().addWidget(clusters_win.getInternalRep())
        self.clusters_win = clusters_win
        clusters_win.getInternalRep().cluster_window = weakref.proxy(self)

        measures_view = QtGui.QWidget()
        measures_view.setObjectName('measures_view')
        main_w.addWidget(measures_view)
        measures_view.setLayout(QtGui.QHBoxLayout())

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
        values = list(self.measurements.keys())
        self.cluster_slider_label = QtGui.QLabel(str(values[0]))
        self.cluster_slider = QtGui.QSlider()
        lay.addWidget(QtGui.QLabel('K:'))
        lay.addWidget(self.cluster_slider_label)
        lay.addWidget(self.cluster_slider)
        self.cluster_slider.setInvertedAppearance(True)
        self.cluster_slider.setInvertedControls(True)
        self.cluster_slider.setRange(values[0], values[-1])
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
            measure_tex.setName(
                list(self.measurements.values())[0].columns[self.viewing_column])
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

        # Hide slider of 2nd window
        # hide both slider and label instead of hiding the parent panel
        # because the panel will be showed automatically by Window3D.Refresh()
        time_slider2 = measures_win.findChild(QtGui.QSlider, 'sliderT')
        time_slider2.hide()
        time_slider2.parentWidget().findChild(QtGui.QLabel).hide()

        self.selected_cluster_boundaries.setMaterial(
            line_width=5., diffuse=[0.6, 0.8, 0., 1.])

    def make_measurements_texture(self):
        col = self.viewing_column
        k = self.cluster_slider.value()
        timestep = k - list(self.measurements.keys())[0]
        for atex, ctex, measure_tex in zip(self.aims_measure_tex,
                                           self.aims_clusters,
                                           self.measure_tex):
            ndata = len(ctex[0].data())
            arr = np.zeros((ndata,), dtype=np.float32)
            if len(ctex) > timestep:
                measurements = self.measurements[k]
                values = measurements[measurements.columns[col]]
                for label in range(len(values)):
                    value = values[label]
                    arr[ctex[timestep].data().arraydata() == label + 1] = value
            atex[0].data().assign(arr)
            measure_tex.setTexture(atex)
            measure_tex.setChanged()
            measure_tex.notifyObservers()


    def print_cluster_info(self, timestep, cluster, mesh=None, ivertex=None):
        color = self.get_cluster_color(timestep, mesh, cluster)
        if color is not None:
            intensity = np.sqrt(np.sum(np.array(color[:3]) ** 2)) / np.sqrt(3.)
            if intensity < 160:
                text_col = ' style="color: white; background-color: ' \
                    '#%02x%02x%02x;"' % color[:3]
            else:
                text_col = ' style="background-color: #%02x%02x%02x;"' \
                    % color[:3]
        else:
            text_col = ''
        info_text = '''<h1%s>Cluster %d info:</h1>
Num of clusters (K): <b>%d</b><br/>
''' \
            % (text_col, cluster, timestep)
        patch = self.get_patch(timestep, mesh, ivertex)
        if patch[0] is not None:
            if patch[1] is not None:
                label = patch[1]['Label']
                color = tuple(int(x * 255.9) for x in patch[1]['RGB'][:3])
                intensity = np.sqrt(np.sum(np.array(color) ** 2)) / np.sqrt(3.)
                if intensity < 135:
                    text_col = ' color: white;'
                else:
                    text_col = ''
                clabel = '<b style="background-color: ' \
                    '#%02x%02x%02x;%s">%s</b>' % (color + (text_col, label))
                info_text += 'In patch (gyrus): %d, label: %s<br/>' \
                    % (patch[0], clabel)
            else:
                info_text += 'In patch (gyrus): %d<br/>' % patch[0]
        self.info.setText(info_text)


    def get_cluster_color(self, timestep, mesh, cluster):
        if mesh is None:
            return None
        mesh_index = self.meshes.index(mesh)
        cluster_tex = self.clusters[mesh_index]
        te = cluster_tex.glTexExtrema()
        pos = (cluster - te.minquant[0]) / (te.maxquant[0] - te.minquant[0])
        pal = cluster_tex.palette()
        colors = pal.colors()
        pos = (pos - pal.min1()) / (pal.max1() - pal.min1())
        if pos < 0:
            pos = 0
        elif pos > 1:
            pos = 1
        color = colors.value(int(pos * (colors.dimX() - 0.0001)))
        color = (color[0], color[1], color[2], color[3])
        return color


    def get_patch(self, timestep, mesh, ivertex):
        if mesh is None or ivertex is None or len(self.seed_gyri) == 0:
            return None, None
        mesh_index = self.meshes.index(mesh)
        seed_gyri = self.aims_seed_gyri[mesh_index]
        patch_label = seed_gyri[0].data()[ivertex]
        hdr = seed_gyri.header()
        if 'GIFTI_labels_table' in hdr:
            label = hdr['GIFTI_labels_table'][patch_label]
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


    def cluster_selected(self, timestep, label, mesh, ivert):
        timestep = self.cluster_slider.value()
        self.print_cluster_info(timestep, label, mesh, ivert)
        self.set_curve_cursor(timestep, label)
        self.select_table_row(timestep, label)
        self.update_cluster_time(timestep, label)
        self.clusters_win.statusBar().showMessage('Cluster: %d' % label)
        self.update_selected_cluster_boundaries(label, mesh)
        self.set_matrix_cursor(label)


    def update_clusters_boundaries(self):
        timestep = self.cluster_slider.value() - list(self.measurements.keys())[0]
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


    def update_selected_cluster_boundaries(self, label, mesh):
        self.clusters_win.removeObjects(self.selected_cluster_boundaries)
        self.measures_win.removeObjects(self.selected_cluster_boundaries)
        if mesh is None or label == 0:
            return
        mesh_index = self.meshes.index(mesh)
        timestep = self.cluster_slider.value() - list(self.measurements.keys())[0]
        cluster_tex0 = self.aims_clusters[mesh_index][timestep]
        cluster_tex = self.aims_clusters[mesh_index].__class__()
        cluster_tex[0].assign(cluster_tex0.data())
        cluster_boundaries = aims.SurfaceManip.meshTextureBoundary(
            mesh.surface(), cluster_tex, label)
        self.selected_cluster_boundaries.setSurface(cluster_boundaries)
        self.selected_cluster_boundaries.setChanged()
        self.selected_cluster_boundaries.notifyObservers()
        self.clusters_win.addObjects(self.selected_cluster_boundaries,
                                     temporary=True, position=0)
        self.measures_win.addObjects(self.selected_cluster_boundaries,
                                     temporary=True, position=0)


    def update_clusters_texture(self):
        timestep = self.cluster_slider.value() - list(self.measurements.keys())[0]
        for clusters, aims_clusters in zip(self.clusters, self.aims_clusters):
            t = min(timestep, len(aims_clusters) - 1)
            texture = aims.TimeTexture('S16')
            texture[0].data().assign(aims_clusters[t].data())
            clusters.setTexture(texture)
            clusters.notifyObservers()


    def time_changed(self, value):
        self.cluster_slider_label.setText(str(value))
        timestep = value
        self.build_table(timestep)
        self.display_curves(timestep)
        self.make_measurements_texture()
        self.update_clusters_texture()
        self.update_clusters_boundaries()
        self.cluster_selected(0, 0, None, 0)
        self.display_matrix()
        #self.update_selected_cluster_boundaries(0, None)


    def display_column(self, col):
        if self.viewing_column != col:
            name = list(self.measurements.values())[0].columns[col]
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
            self.curves_fig.subplots_adjust(bottom=0.2, right=0.8)
        else:
            # use existing axes
            axes = self.curves_fig.axes[0]
            axes.clear()
        axes.set_xlabel('cluster')
        data = self.measurements[timestep]
        cols = data.columns[self.curves_columns]
        if len(cols) != 0:
            axes.set_xticks(list(range(data.shape[0] + 1)))
            axes.set_xlim((0.8, data.shape[0] + 0.2))
            axes.plot(list(range(1, data.shape[0] + 1)), data[cols], 'o-')
            axes.legend(list(cols), loc=(1.02, 0.2), fontsize='x-small',
                        labelspacing=0.5)
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
                      lim + lim[::-1] + (lim[0], ), 'r', linewidth=2)
            self.curves_fig.canvas.draw()


    def select_curves_columns(self):
        dialog = QtGui.QDialog()
        dialog.setModal(True)
        lay = QtGui.QVBoxLayout()
        dialog.setLayout(lay)
        lwid = QtGui.QListWidget()
        lay.addWidget(lwid)
        lwid.addItems(list(self.measurements.values())[0].columns)
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
            timestep = self.cluster_slider.value()
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


    def build_matrix_by_cluster(self, k):
        if self.matrix is None:
            return None
        if k in self.clusters_matrix:
            return self.clusters_matrix[k]
        timestep = k - list(self.measurements.keys())[0]
        matrix = self.matrix
        a_mat = np.asarray(matrix)
        a_mat = a_mat.reshape(a_mat.shape[:2])
        clusters_arr = self.aims_clusters[0][timestep].arraydata() ## WARNING use tex 0
        new_matrix = matrixtools.compute_mclusters_by_nbasins_matrix(
            a_mat, clusters_arr, mode='mean')
        self.clusters_matrix[k] = new_matrix
        return new_matrix


    def display_matrix(self):
        k = self.cluster_slider.value()
        matrix = self.build_matrix_by_cluster(k)
        if len(self.matrix_fig.axes) == 0:
            # no axes yet. Create them
            bgcolor = self.palette().color(QtGui.QPalette.Base)
            axes = self.matrix_fig.add_subplot(
                111, axisbg=str(bgcolor.name()))
            self.matrix_fig.subplots_adjust(bottom=0.2)
        else:
            axes = self.matrix_fig.axes[0]
            axes.clear()
        axes.set_xlabel('basins')
        axes.set_ylabel('clusters')
        axes.set_yticks(list(range(matrix.shape[0] + 1)))
        axes.set_yticklabels([str(x) for x in range(1, matrix.shape[0] + 1)])
        axes.imshow(matrix, aspect='auto', interpolation='none')
        self.matrix_fig.canvas.draw()


    def set_matrix_cursor(self, label):
        if len(self.curves_fig.axes) != 0:
            axes = self.matrix_fig.axes[0]
            if len(axes.lines) != 0:
                # erase previous cursor
                axes.lines.remove(axes.lines[-1])
            if label == 0:
                return
            lim = axes.get_xlim()
            xlim = (lim[0] + 0.01, lim[1] - 0.01)
            ylim = (label - 1.49, label - 0.51)
            axes.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]],
                      [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], 'r',
                      linewidth=2)
            axes.set_xlim(lim)
            self.matrix_fig.canvas.draw()


def load_measurements(measurements_filename):
    measurements = {}
    with open(measurements_filename) as f:
        title = f.readline().strip()[1:-1]
        columns = eval(title)
        lnum = 1
        #k0 = 0
        for l in f:
            lnum += 1
            t = l.find(',')
            k = int(l[:t].strip())
            #if k0 == 0:
                #k0 = k
            table = np.zeros((k, len(columns)))
            meas = l[t+1:].strip()
            for i in range(k):
                if meas[0] != '"':
                    raise SyntaxError('missing double quote, line: %d' % lnum)
                t = meas[1:].find('"')
                m = eval(meas[1:t+1])
                table[i,:len(m)] = m
                meas = meas[t+2:].strip()
                if meas:
                    if meas[0] != ',':
                        raise SyntaxError('missing coma, line: %d' % lnum)
                    meas = meas[1:]
            measurements[k] = pandas.DataFrame(table, columns=columns)
    #print('measurements:', measurements)
    return measurements


def load_clusters_inspector_files(mesh_filenames, clusters_filenames,
                                  measurements_filenames,
                                  seed_gyri_filenames=[],
                                  matrix_filename=None):
    if len(mesh_filenames) != len(clusters_filenames):
        raise ValueError('meshes and clusters numbers do not match')
    if len(seed_gyri_filenames) != 0 \
            and len(seed_gyri_filenames) != len(mesh_filenames):
        raise ValueError('meshes and seed_gyri numbers do not match')

    meshes = []
    clusters = []
    seed_gyri = []
    matrix = None
    measurements = {}
    for mesh_filename in mesh_filenames:
        meshes.append(aims.read(mesh_filename))
    for clusters_filename in clusters_filenames:
        clusters.append(aims.read(clusters_filename))
    if measurements_filenames:
        for measurements_filename in measurements_filenames:
            measurements.update(load_measurements(measurements_filename))
    for seed_gyri_filename in seed_gyri_filenames:
        seed_gyri.append(aims.read(seed_gyri_filename))

    if matrix_filename is not None:
        matrix = aims.read(matrix_filename)

    return meshes, clusters, measurements, seed_gyri, matrix


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

    use_ex_num = 2

    if use_ex_num == 0:
        meshes, clusters, measurements, seed_gyri, matrix \
            = load_clusters_inspector_files(
                #['/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.3.0/group_analysis/01to40/average_brain/averagebrain.white.mesh'],
                ['/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.1.0/001/surf/bh.r.aims.white.inflated.gii'],
                ['/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to40/connectivity_clustering/avg/fs01to40/lh.supramarginal/smooth3.0/avgSubject/01to40_avg_fs01to40_lh.supramarginal_avgSubject_clusteringTime.gii'],
                None,
                #['/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.3.0/group_analysis/01to40/average_brain/bh.annot.averagebrain.gii'])
                ['/tmp/bh.r.aparc.annot.gii'],
                '/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to40/connectivity_clustering/avg/fs01to40/lh.supramarginal/smooth3.0/01to40_avg_fs01to40_lh.supramarginal_matrix.ima')
    elif use_ex_num == 1:
        meshes, clusters, measurements, seed_gyri, matrix \
            = load_clusters_inspector_files(
                ['/neurospin/population/HCP/S500-1/100307/T1w/fsaverage_LR32k/100307.L.inflated.32k_fs_LR.surf.gii', '/neurospin/population/HCP/S500-1/100307/T1w/fsaverage_LR32k/100307.R.inflated.32k_fs_LR.surf.gii'],
                ['/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to40/connectivity_clustering/avg/fs01to40/lh.supramarginal/smooth3.0/avgSubject/01to40_avg_fs01to40_lh.supramarginal_avgSubject_clusteringTime.gii', '/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to40/connectivity_clustering/avg/fs01to40/lh.supramarginal/smooth3.0/avgSubject/01to40_avg_fs01to40_lh.supramarginal_avgSubject_clusteringTime.gii'],
                None,
                ['/neurospin/population/HCP/S500-1/100307/MNINonLinear/fsaverage_LR32k/100307.L.aparc.32k_fs_LR.label.gii', '/neurospin/population/HCP/S500-1/100307/MNINonLinear/fsaverage_LR32k/100307.R.aparc.32k_fs_LR.label.gii'])
    elif use_ex_num == 2:
        meshes, clusters, measurements, seed_gyri, matrix \
            = load_clusters_inspector_files(
                ['/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.3.0/001/surf/bh.r.aims.white.inflated_sequence.gii'],
                ['/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to39/connectivity_clustering/avg/fs01to39/lh.postcentral/smooth3.0/avgSubject/01to39_avg_fs01to39_lh.postcentral_avgSubject_clusteringTime.gii'],
                None,
                ['/neurospin/archi-public/DataBaseArchi/FreeSurfer/fs_archi_v5.3.0/group_analysis/01to39/average_brain/bh.annot.averagebrain.gii'],
                '/neurospin/archi-public/Users/lefranc/archi/bv_archi/proba27/subjects/group_analysis/01to39/connectivity_clustering/avg/fs01to39/lh.postcentral/smooth3.0/01to39_avg_fs01to39_lh.postcentral_matrix.ima',
            )
    else:
        print('wrong use_ex_num value:', use_ex_num)
        raise ValueError('aborting')
    # temp
    measurements = dict(
        (i, pandas.DataFrame(np.random.ranf((i + 2, 4)),
                             columns=('size', 'homogeneity', 'conn_density',
                                      'other')))
        for i in range(max([len(clusters_tex) for clusters_tex in clusters])))
    cw = ClustersInspectorWidget(
        meshes, clusters, measurements=measurements, seed_gyri=seed_gyri,
        matrix=matrix)
    cw.show()

    if run_event_loop:
        QtGui.QApplication.instance().exec_()


