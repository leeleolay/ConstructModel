#!/bin/bash

# 声明包含所有文件夹的数组
folders_old=(
    "230814_1.1.1gamma_inter"
    "230907_1.1.2gamma_p_gamma_inter"
    "230908_1.2.1a_inter"
    "230908_1.2.2a_p_a_inter"
    "230909_2.1temp_mu_sigma"
    "230909_2.2temp_sigma"
    "230909_3.1Force"
    "230910_3.2.1velocitymagnitude"
    "230910_3.2.2velocitydirection"
    "230910_4.1cellgap"
)

folders=(
	"231016_1.1.1gamma_inter"
	"231016_1.1.2gamma_p_gamma_inter"
)

# 循环遍历每个文件夹
for folder in "${folders[@]}"; do
    # 构建完整的文件夹路径
    folder_path="/ncsfs02/lxw/taskCurrent/${folder}"

    # 获取文件夹中以 "lammpstrj." 为前缀的数据文件名列表
    data_files=$(find "${folder_path}" -maxdepth 1 -type f -name "lammpstrj.*")

    # 打印文件名列表
    echo "Data Files in ${folder}:"
    echo "${data_files}"

    # 循环遍历每个数据文件
    for data_file in ${data_files}; do
        # 构建 TransportTime.py 命令
        transport_time_cmd="python3 /ncsfs02/lxw/opt/Toolkit4MesoCellModel/PostProcess/TransportTime.py --file_name ${data_file} --atom_type 3 --min_timestep 0 --max_timestep 500000 --particle_ids 97409 97410 97411 97412 97413 97414 97415 97416 97417 97418 97419 97420 97421 97422 97423 97424 97425 97426 97427 97428 97429 97430 97431 97432 97433 97434 97435 97436 97437 97438 97439 97440 97441 97442 97443 97444 97445 97446 97447 97448 97449 97450 97451 97452 97453 97454 97455 97456 97457 97458 97459 97460 97461 97462 97463 97464 97465 97466 97467 97468 97469 97470 97471 97472 97473 97474 97475 97476 97477 97478 97479 97480 97481 97482 97483 97484 97485 97486 97487 97488 97489 97490 97491 97492 97493 97494 97495 97496 97497 97498 97499 97500 97501 97502 97503 97504 97505 97506 97507 97508 97509 97510 97511 97512 97513 97514 97515 97516 97517 97518 97519 97520 97521 97522 97523 97524 97525 97526 97527 97528 97529 97530 97531 97532 97533 97534 97535 97536 97537 97538 97539 97540 97541 97542 97543 97544 97545 97546 97547 97548 97549 --start_region -42 -42 25.0 42 42 35.0 --end_region -42 -42 -35.0 42 42 -25.0 --output ${folder_path}/ "

        # 打印 TransportTime.py 命令
        echo "Running TransportTime.py for file: ${data_file}"
        echo "Command: ${transport_time_cmd}"

        # 执行 TransportTime.py 命令
        eval ${transport_time_cmd}

        echo "TransportTime.py completed for file: ${data_file}"

        # 构建 MSD.py 命令
        msd_cmd="python3 /ncsfs02/lxw/opt/Toolkit4MesoCellModel/PostProcess/MSD.py ${data_file} 3 0 11300 ${folder_path}/new/ --num_processes 14"

        # 打印 MSD.py 命令
        echo "Running MSD.py for file: ${data_file}"
        echo "Command: ${msd_cmd}"

        # 执行 MSD.py 命令
        #eval ${msd_cmd}

        echo "MSD.py completed for file: ${data_file}"
    done
done
