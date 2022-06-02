from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing, QgsFeatureSink, QgsProcessingException, QgsProcessingAlgorithm, QgsProcessingParameterRasterLayer, QgsProcessingParameterBoolean, QgsProcessingParameterRasterDestination)
from qgis import processing

from topopy import DEM, Flow

class FlowDirection(QgsProcessingAlgorithm):
    INPUT_DEM = 'INPUT_DEM' #Input DEM grid
    OUTPUT_FD = 'OUTPUT_FD' #Output flow direction raster
    VERBOSE = 'VERBOSE'
    AUXTOPO = 'AUXTOPO'
    FILLED = 'FILLED'
        
    def tr(self, string):
        return QCoreApplication.translate('flowdirection', string)

    def createInstance(self):
        return type(self)()

    def name(self):
        return 'flowdirection'

    def displayName(self):
        return self.tr('Flow Direction')

    def group(self):
        return self.tr('topopy')

    def groupId(self):
        return 'topopy'

    def shortHelpString(self):
        texto = "Takes a DEM as inpuut and outputs a flow direction raster produced by an object of the Flow class."
        return texto

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DEM,  self.tr("Input DEM")))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_FD, "Output Flow Direction raster", None, False))
        self.addParameter(QgsProcessingParameterBoolean(self.AUXTOPO, "Use auxiliar topography (ignored if Filled DEM is set to True)", False))
        self.addParameter(QgsProcessingParameterBoolean(self.FILLED, "Filled DEM", False))
        self.addParameter(QgsProcessingParameterBoolean(self.VERBOSE, "Show messages", False))
        
    def processAlgorithm(self, parameters, context, feedback):
        dem_raster = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        fd_raster = self.parameterAsOutputLayer(parameters, self.OUTPUT_FD, context)
        verbose = self.parameterAsBool(parameters, self.VERBOSE, context)
        filled = self.parameterAsBool(parameters, self.FILLED, context)
        auxtopo = self.parameterAsBool(parameters, self.AUXTOPO, context)
        
        dem = DEM(dem_raster.source()) #Load DEM
        fd = Flow(dem, auxtopo=auxtopo, filled=filled, verbose=verbose, verb_func=feedback.setProgressText) #Create Flow object
        fd.save(fd_raster) #Save Flow object onto disc
        results = {self.OUTPUT_FD : fd_raster}
        return results
