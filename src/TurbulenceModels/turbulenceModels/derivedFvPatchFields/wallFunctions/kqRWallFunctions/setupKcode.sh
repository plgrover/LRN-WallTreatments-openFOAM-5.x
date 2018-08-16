#!/bin/bash

clear

cp -r --parents src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/kqRWallFunctions/kLowReWallFunction $WM_PROJECT_USER_DIR 

mv kLowReWallFunction kLowReColebrookWallFunction
cd kLowReColebrookWallFunction
rm *.dep

mv kLowReWallFunctionFvPatchScalarField.C kLowReColebrookWallFunctionFvPatchScalarField.C
mv kLowReWallFunctionFvPatchScalarField.H kLowReColebrookWallFunctionFvPatchScalarField.H

sed -i s/kLowReWallFunction/kLowReColebrookWallFunction/g kLowReColebrookWallFunctionFvPatchScalarField.C
sed -i s/kLowReWallFunction/kLowReColebrookWallFunction/g kLowReColebrookWallFunctionFvPatchScalarField.H

echo "done"

