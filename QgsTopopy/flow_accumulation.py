from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing, QgsFeatureSink, QgsProcessingException, QgsProcessingAlgorithm, QgsProcessingParameterRasterLayer, QgsProcessingParameterBoolean, QgsProcessingParameterRasterDestination)
from qgis import processing

from topopy import DEM, Flow

class FlowAccumulation(QgsProcessingAlgorithm):
    INPUT_DEM = 'INPUT_DEM' #Input DEM grid
    OUTPUT_FAC = 'OUTPUT_FAC' #Output flow accumulation raster
    VERBOSE = 'VERBOSE'
    AUXTOPO = 'AUXTOPO'
    FILLED = 'FILLED'
        
    def tr(self, string):
        return QCoreApplication.translate('flowaccumulation', string)

    def createInstance(self):
        return type(self)()

    def name(self):
        return 'flowaccumulation'

    def displayName(self):
        return self.tr('Flow Accumulation')

    def group(self):
        return self.tr('topopy')

    def groupId(self):
        return 'topopy'

    def shortHelpString(self):
        texto = "Takes a DEM as inpuut and outputs a flow accumulation raster using the Flow.get_flow_accumulation() method. Keeps default parameters (mantains NoData values, outputs a topopy.Grid object)."
        return texto

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DEM,  self.tr("Input DEM")))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_FAC, "Output Flow Accumulation raster", None, False))
        self.addParameter(QgsProcessingParameterBoolean(self.AUXTOPO, "Use auxiliar topography (ignored if Filled DEM is set to True)", False))
        self.addParameter(QgsProcessingParameterBoolean(self.FILLED, "Filled DEM", False))
        self.addParameter(QgsProcessingParameterBoolean(self.VERBOSE, "Show messages", False))
        
    def processAlgorithm(self, parameters, context, feedback):
        dem_raster = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        fac_raster = self.parameterAsOutputLayer(parameters, self.OUTPUT_FAC, context)
        verbose = self.parameterAsBool(parameters, self.VERBOSE, context)
        filled = self.parameterAsBool(parameters, self.FILLED, context)
        auxtopo = self.parameterAsBool(parameters, self.AUXTOPO, context)
        
        dem = DEM(dem_raster.source()) #Load DEM
        fd = Flow(dem, auxtopo=auxtopo, filled=filled, verbose=verbose, verb_func=feedback.setProgressText) #Create Flow object
        fac = fd.get_flow_accumulation() #Get flow accumulation raster
        fac.save(fac_raster) #Save flow accumulation raster onto disc
        results = {self.OUTPUT_FAC : fac_raster}
        return results
