#/bin/bash
export http_proxy=
export https_proxy=
curl –v  "localhost:8000/api/workflows/v1" -F workflowSource=@workflows/Optimized_20k_Single_Sample.wdl -F workflowInputs=@workflows/Optimized_20k_Single_Sample_inputs.json > 20k_submission_response.txt
cat submission_response.txt | grep -o -E "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}" > 20k_WF_ID.txt

