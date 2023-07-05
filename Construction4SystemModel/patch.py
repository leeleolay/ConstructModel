import os

# 根目录，这里替换成你的根目录路径
root_dir = os.getcwd()

# Python脚本路径
python_script = os.path.join(root_dir, "transform_to_lammps.py")

# 创建一个列表，存放所有的子目录
all_dirs = []

# 遍历根目录下的所有子目录，添加到列表中
for subdir, dirs, files in os.walk(root_dir):
    for dir in dirs:
        # 拼接完整的子目录路径
        all_dirs.append(os.path.join(subdir, dir))

# 遍历所有的子目录，执行脚本
for idx, model_dir in enumerate(all_dirs):
    # 保存当前工作目录
    current_dir = os.getcwd()

    # 进入目标目录
    os.chdir(model_dir)

    # 打印进度信息
    print(f"Processing {idx+1} of {len(all_dirs)}: {model_dir}")

    # 执行你的Python脚本，当一个子进程结束后，再开始下一个子进程
    with open(python_script) as file:
        exec(file.read())

    # 返回到原来的工作目录
    os.chdir(current_dir)


