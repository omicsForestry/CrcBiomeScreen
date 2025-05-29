# CrcBiomeScreen

![Pipeline](https://private-user-images.githubusercontent.com/82665896/448789645-51279bdb-eeff-45bc-bc46-2b8a72c43a5a.jpg?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NDg1MTQ1MTAsIm5iZiI6MTc0ODUxNDIxMCwicGF0aCI6Ii84MjY2NTg5Ni80NDg3ODk2NDUtNTEyNzliZGItZWVmZi00NWJjLWJjNDYtMmI4YTcyYzQzYTVhLmpwZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA1MjklMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwNTI5VDEwMjMzMFomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTAxYTU4YmU2OGMxMGYwYTNkNjRmNDhkMjk0NzViOGVhZTA5M2QxOGMwN2JiOTk5ZDBhYTlkYTQ0NmIxNTQ3YzYmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.YMfd4r5CyLUizv0jCyzbbVApW2oN8Z2cQDeJOSqEWCE)

## Description
 * Vignette.R : The toy sample for this workflow.
 * Environment.R : For setting up the environment for this workflow.
 * R folder : Multiple functions for processing and analyzing.
 * CrcBiomeScreen.png

## Memos...
 * Need to pay attention to the <mark>'class_weights'</mark> in differnt datasets.
   * How to define the imbalance ratio better to improve applicability...[Part of it]
 * Need to add the <mark>attributes</mark> in the object![Maybe no need to do this part...]

## Plan
 * Need to do the test for the whole pipeline,especially the PredictValidation part.✅
 * Add the description and annotation in Vignette and function.➡️
 * Record all operation history ➡️
 * In the future, more functions(modules) will put in this workflow.
   * CompareModels()
   * SelectImportanceFeatures()
   * SHAP values...
