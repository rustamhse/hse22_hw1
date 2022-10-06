# Биоинформатика 2022 ДЗ-1
###### Сборка генома бактерии, выделенной из воды с нефтью, на основании парно-концевых (paired-end, PE) и чтений типа mate-pairs (MP). 
###### Зиёев Рустам Группа 3

## Работа с Putty

1. Создание символических ссылок
``` 
$ ls /usr/share/data-minor-bioinf/assembly/* | xargs -tI{} ln -s {} 
```

2. Выбор случайных чтений с помощью SeqTK (Random seed = 928)
``` $seqtk sample -s 928 oil_R1.fastq 5000000 > sub1.fastq
$ seqtk sample -s 928 oil_R2.fastq 5000000 > sub2.fastq
$ seqtk sample -s 928 oilMP_S4_L001_R1_001.fastq 1500000 > matep1.fastq
$ seqtk sample -s 928 oilMP_S4_L001_R2_001.fastq 1500000 > matep2.fastq
```

3. Оценка качества исходных чтений с помощью FastQC и создание отчёта с помощью MultiQC
```
$ mkdir fastqc
$ ls sub* matep* | xargs -tI{} fastqc -o fastqc {}
$ mkdir multiqc
$ multiqc -o multiqc fastqc
```

4. Урезание размера исходных чтений через Platanus, оценка и создание отчёта MultiQC:
```
$ platanus_trim sub*
$ platanus_internal_trim matep*
$ mkdir fastqc_trimmed
$ ls sub* matep*| xargs -tI{} fastqc -o fastqc_trimmed {}
$ mkdir multiqc_trimmed
$ multiqc -o multiqc_trimmed fastqc_trimmed
```

5. Сбор контиг и скаффолдов

```
$ platanus assemble -o Poil -f sub1.fastq.trimmed sub2.fastq.trimmed
$ platanus scaffold -o Poil -c Poil_contig.fa -IP1 sub1.fastq.trimmed sub2.fastq.trimmed -OP2 matep1.fastq.int_trimmed matep2.fastq.int_trimmed
$ platanus gap_close -o Poil -c Poil_scaffold.fa -IP1 sub1.fastq.trimmed sub2.fastq.trimmed -OP2 matep1.fastq.int_trimmed matep2.fastq.int_trimmed
```

## Отчёты

Для неурезанных чтений

![image](https://user-images.githubusercontent.com/74643940/194400526-fc75a017-efa3-4007-914a-edea8c7b1178.png)

![image](https://user-images.githubusercontent.com/74643940/194401775-678664a6-9327-4b8c-8cb0-b69a0ae2c5a5.png)


Для урезанных чтений

![image](https://user-images.githubusercontent.com/74643940/194401439-54e4408b-5983-49aa-a57a-392e29e21c59.png)

![image](https://user-images.githubusercontent.com/74643940/194401841-4c6da966-d660-4c52-bade-fb8d754a1e8b.png)

## Colab
###### https://colab.research.google.com/drive/1KD7XxoKcSD4SAvnd5Pcl5zljIr9w-ajq?usp=sharing

**Инициализация функции для вывода информации**

```
import re

def get_info(f, output = True):
    l = []
    num = 0 
    max_len = 0
    length = 0
    t_l = 0
    s = 0
    max_s = ''
    curr_s = ''
    
    for line in f:
        if (line[0] == '>'):
            if num != 0:
                l.append(length)
            num += 1
            if length >= max_len:
                max_len = length
                max_s = curr_s
            curr_s = ''
            length = 0
        else:
            curr_s += line.strip()
            length += len(line.strip())
            t_l += len(line.strip())
     
    l.sort(reverse = True) 
    for i in l:
        s += i
        if s >= t_l / 2:
            if output == True:
                print(f'Results:\n\
                Overall quantity: {num},\n\
                Overall length: {t_l},\n\
                The longest section: {max_len},\n\
                N50 mean: {i}\n')
            break
    return max_s```
```

**Контиги**

Platanus assemble

```
max_cont = get_info(open('Poil_contig.fa', 'r'))
```

>Results:
>Overall quantity: 612,
>Overall length: 3925328,
>The longest section: 179307,
>N50 mean: 47611

**Скаффолды**

Platanus scaffold

```
max_scaf = get_info(open('Poil_scaffold.fa', 'r'))
```

>Results:
>Overall quantity: 70,
>Overall length: 3876568,
>The longest section: 3838373,
>N50 mean: 3838373


```
print(f'Overall gaps length: {max_scaf.count("N")}')
max_scaf = re.sub(r'N{2,}', 'N', max_scaf)
print(f'Gaps quantity: {max_scaf.count("N")}')
```

>Overall gaps length: 6346
>Gaps quantity: 65


```
m_s = get_info(open('Poil_gapClosed.fa', 'r'), False)
m_s_sub = re.sub(r'N{2,}', 'N', m_s)
print(f'Overall sub gaps length: {m_s.count("N")}')
print(f'Sub gaps quantity: {m_s_sub.count("N")}')
```

>Overall sub gaps length: 1674
>Sub gaps quantity: 9

