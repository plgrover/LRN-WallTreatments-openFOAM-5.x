clear


mv omegaWallFunctionFvPatchScalarField.C omegaTestWallFunctionFvPatchScalarField.C
mv omegaWallFunctionFvPatchScalarField.H omegaTestWallFunctionFvPatchScalarField.H

sed -i s/omegaWallFunction/omegaColebrookWallFunction/g omegaTestWallFunctionFvPatchScalarField.C
sed -i s/omegaWallFunction/omegaColebrookWallFunction/g omegaTestWallFunctionFvPatchScalarField.H

echo "done"
