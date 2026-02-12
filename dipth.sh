docker run -it --rm \
    -v /home/user1/Desktop/programs/wdl/IsraBact/data:/data \
    bioinfomoh/diphtoscan:latest \
    /bin/bash -c "source activate diphtoscan && diphtoscan \
        -a /home/user1/Desktop/programs/wdl/IsraBact/data/CB521_filtered_contigs.fasta /home/user1/Desktop/programs/wdl/IsraBact/data/CB522_filtered_contigs.fasta \
        -st -t -res_vir -plus --integron \
        -o /data/results --threads 2 "
