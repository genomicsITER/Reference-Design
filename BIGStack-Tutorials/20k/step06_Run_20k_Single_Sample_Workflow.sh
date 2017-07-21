#/bin/bash
export http_proxy=
export https_proxy=
curl -XPOST -u USERNAME:PASSWORD http://localhost:8080/api/v1/workflowcollections -F workflowSource=@workflows/Optimized_20k_Single_Sample.wdl -F workflowInputs=@workflows/Optimized_20k_Single_Sample_inputs.json > 20k_submission_response.txt
cat submission_response.txt | grep -o -E "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}" > 20k_WF_ID.txt

