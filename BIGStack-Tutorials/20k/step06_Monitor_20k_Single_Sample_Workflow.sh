#/bin/bash
export http_proxy=
export https_proxy=
curl –v  "localhost:8000/api/workflows/v1/$(cat 20k_WF_ID.txt)/status

