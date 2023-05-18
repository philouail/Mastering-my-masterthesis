setImageType('FLUORESCENCE');


def pathModel = "D:\\Users\\zeiss\\Desktop\\Phili\\scripts supporting open-source tools for PhenoCycler image analysis-2\\stardist_cell_seg_model.pb"
import qupath.ext.stardist.StarDist2D

def stardist = StarDist2D.builder(pathModel)
        .threshold(0.5)              // Probability (detection) threshold
        .channels(6)            // Select detection channel
        .normalizePercentiles(1, 99) // Percentile normalization
        .pixelSize(0.25)              // Resolution for detection
        .cellExpansion(5.0)          // Approximate cells based upon nucleus expansion
        .cellConstrainScale(1.5)     // Constrain cell expansion using nucleus size
        .measureShape()              // Add shape measurements
        .measureIntensity()          // Add cell measurements (in all compartments)
        .includeProbability(true)    // Add probability as a measurement (enables later filtering)
        .build()

// Run detection for the selected objects
def imageData = getCurrentImageData()
def pathObjects = getSelectedObjects()
createSelectAllObject(true)
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
    return
}
stardist.detectObjects(imageData, pathObjects)
println 'Done segmentation !'

createAnnotationsFromPixelClassifier("PDL1_4", 0.0, 0.0)
selectCells();
println 'Done creation of PDL1 object !'


addPixelClassifierMeasurements("PDL1_4", "PDL1_4")
resetSelection();
println 'Done measurement of PDL1 !'


createAnnotationsFromPixelClassifier("Tumor_stroma", 0.0, 0.0)
selectCells();
println 'Done Tumor and stroma obejct creation !'

addPixelClassifierMeasurements("Tumor_stroma", "Tumor_stroma")
println 'Done Tumor and Stroma measurement !'

selectObjectsByClassification("Stroma");
createAnnotationsFromPixelClassifier("Fibroblast", 0.0, 0.0)
selectCells();
println 'Done Fibroblast object within stroma !'

addPixelClassifierMeasurements("Fibroblast", "Fibroblast")
runObjectClassifier("Immune_cells_2");
println 'Done fibroblast measurement, end of script !'

//DO UNKNOWN ANALYSIS SEPARATELY
