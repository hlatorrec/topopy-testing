from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing, QgsProcessingException, QgsFeatureSink, QgsProcessingAlgorithm,
                       QgsProcessingParameterRasterLayer, QgsProcessingParameterFileDestination,
                       QgsProcessingParameterNumber, QgsProcessingParameterBoolean)
from qgis import processing

from topopy import Flow, Network


class GetNetwork(QgsProcessingAlgorithm):
    INPUT_FD = 'INPUT_FD'
    OUTPUT_NETW = 'OUTPUT_NETW'
    THRESHOLD = 'THRESHOLD'
    THETAREF = 'THETAREF'
    GRADIENTS = 'GRADIENTS'

    def tr(self, string):
        return QCoreApplication.translate('getnetwork', string)

    def createInstance(self):
        return type(self)()

    def name(self):
        return 'getnetwork'

    def displayName(self):
        return self.tr('Get Network')

    def group(self):
        return self.tr('topopy')

    def groupId(self):
        return 'topopy'

    def shortHelpString(self):
        return self.tr("")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_FD, self.tr('Input Flow Direction raster')))
        self.addParameter(QgsProcessingParameterFileDestination(self.OUTPUT_NETW, self.tr('Path to save Network object'), fileFilter = '.dat'))
        self.addParameter(QgsProcessingParameterNumber(self.THRESHOLD, self.tr('Number of cells to initiate a channel'), defaultValue = 0))
        self.addParameter(QgsProcessingParameterNumber(self.THETAREF, self.tr('m/n coeficient to calculate chi values in each channel cell'), QgsProcessingParameterNumber.Double, defaultValue = 0.45))
        self.addParameter(QgsProcessingParameterBoolean(self.GRADIENTS, self.tr('Calculate gradients'), False, True))
        
    def processAlgorithm(self, parameters, context, feedback):
        fd_raster = self.parameterAsRasterLayer(parameters, self.INPUT_FD, context)
        netw_file = self.parameterAsFileOutput(parameters, self.OUTPUT_NETW, context)
        threshold = self.parameterAsInt(parameters, self.THRESHOLD, context)
        thetaref = self.parameterAsDouble(parameters, self.THETAREF, context)
        gradients = self.parameterAsBool(parameters, self.GRADIENTS, context)

        fd = Flow(fd_raster.source()) #Load Flow Direction raster
        netw = Network(fd, threshold = threshold, thetaref = thetaref, gradients = gradients)
        netw.save(netw_file) #Save Network object onto disc

        return {}
