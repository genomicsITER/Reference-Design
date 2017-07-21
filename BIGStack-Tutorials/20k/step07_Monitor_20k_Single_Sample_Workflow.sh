#/bin/bash
export http_proxy=
export https_proxy=
curl -u -XGET USERNAME:PASSWORD 'http://localhost:8080/api/v1/workflows/$(cat 20k_WF_ID.txt)/metadata'

