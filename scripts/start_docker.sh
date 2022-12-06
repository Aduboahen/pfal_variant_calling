# start docker container

docker run \
	-i \
	-t \
	--rm --name pfal_res \
	--mount type=bind,source=${PWD}/input,target=/tmp/input \
        --mount type=bind,source=${PWD}/output,target=/tmp/output \
        --mount type=bind,source=${PWD}/scripts,target=/tmp/scripts \
       	mugicrow/pfal_variant_call:latest

       
