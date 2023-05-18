//Change for the classifier you have created for the specific unknown

createAnnotationsFromPixelClassifier("SLAMF7", 0.0, 0.0)
selectCells();
println 'Done creation of Test Marker object !'


addPixelClassifierMeasurements("SLAMF7", "SLAMF7")
resetSelection();
println 'Done measurement of Test Marker !'