import matplotlib.pyplot as plt

def plot_interpolation_results():
    data = {
        'interpolation_results.txt': 'Interpolation Results',
        'chebyshev_interpolation_results.txt': 'Chebyshev Interpolation Results',
        'heart_curve_10.txt': 'Heart Shape Approximation (10 Segments)',
        'heart_curve_40.txt': 'Heart Shape Approximation (40 Segments)',
        'heart_curve_160.txt': 'Heart Shape Approximation (160 Segments)'
    }

    for filename, title in data.items():
        # 读取数据文件
        x, y = [], []
        with open(filename, 'r') as file:
            for line in file:
                # 移除前后空格并跳过空行和非数值行
                line = line.strip()
                if not line or not line[0].isdigit() and line[0] not in ('-', '+'):
                    continue  # 跳过非数值行

                try:
                    if filename in ['interpolation_results.txt', 'chebyshev_interpolation_results.txt']:
                        _, x_val, y_val = map(float, line.split())
                    else:
                        x_val, y_val = map(float, line.split()[:2])
                    x.append(x_val)
                    y.append(y_val)
                except ValueError:
                    print(f"Skipping invalid line in {filename}: {line}")
                    continue


        # 绘制图表
        plt.figure()
        plt.plot(x, y, label=title, marker='o', markersize=3)
        plt.title(title)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.legend()
        plt.grid(True)
        
        # 保存为图片
        image_filename = filename.replace('.txt', '.png')
        plt.savefig(image_filename)
        plt.close()
        print(f"Saved plot to {image_filename}")

plot_interpolation_results()