
```bash
# bash脚本

#!/bin/bash

# 定义原始蛋白序列
PROTEIN_SEQ="MSTRSVSSSSYRRMFGGPGTASRPSSSRSYVTTSTRTYSLGSALRPSTSRSLYASSPGGVYATRSSAVRL"

# 定义输出文件夹
OUTPUT_DIR="blast_results"
mkdir -p $OUTPUT_DIR

# 随机打乱序列并生成10个新序列
for i in {1..10}; do
    # 使用 `fold` 将序列分割成单个字符，然后用 `shuf` 打乱顺序
    SHUFFLED_SEQ=$(echo $PROTEIN_SEQ | fold -w1 | shuf | tr -d '\n')
    echo ">Sequence_$i" > $OUTPUT_DIR/sequence_$i.fasta
    echo $SHUFFLED_SEQ >> $OUTPUT_DIR/sequence_$i.fasta
done

# 使用 blastp 进行两两比对
for i in {1..10}; do
    for j in {1..10}; do
        if [ $i -ne $j ]; then
            # 比对序列并输出结果
            blastp -query $OUTPUT_DIR/sequence_$i.fasta -subject $OUTPUT_DIR/sequence_$j.fasta -outfmt 6 > $OUTPUT_DIR/blast_result_${i}_${j}.txt
            echo "Blast comparison between Sequence_$i and Sequence_$j saved to $OUTPUT_DIR/blast_result_${i}_${j}.txt"
        fi
    done
done

echo "All blast comparisons completed. Results saved in $OUTPUT_DIR."
```

```bash

# 输出结果文件
test@bioinfo_docker:~/blast/blast_results$ cat blast_result_*.txt
Sequence_10	Sequence_3	35.000	40	22	1	14	53	12	47	1.4	11.5
Sequence_10	Sequence_4	55.556	9	4	0	31	39	15	3.0	10.8
Sequence_10	Sequence_5	41.667	24	14	0	33	56	45	68	0.27	13.5
Sequence_10	Sequence_5	50.000	8	4	0	9	16	58	65	1.8	11.2
Sequence_10	Sequence_5	50.000	10	5	0	30	39	59	68	5.9	10.0
Sequence_10	Sequence_8	48.148	27	12	1	20	44	29	0.022	16.5
Sequence_10	Sequence_9	50.000	18	9	0	39	56	26	0.16	14.2
Sequence_10	Sequence_9	45.455	11	6	0	23	33	59	69	2.5	10.8
Sequence_1	Sequence_2	42.857	28	13	2	36	63	37	61	0.093	14.6
Sequence_1	Sequence_2	42.857	21	9	1	4	24	38	55	0.77	12.3
Sequence_1	Sequence_3	53.846	13	6	0	55	67	21	33	1.3	11.5
Sequence_1	Sequence_3	50.000	6	3	0	38	43	17	22	8.2	 9.6
Sequence_1	Sequence_3	75.000	4	1	0	60	63	57	60	8.8	 9.2
Sequence_1	Sequence_6	50.000	10	5	0	9	18	44	53	0.43	13.1
Sequence_1	Sequence_6	83.333	6	1	0	62	67	21	26	0.78	12.3
Sequence_1	Sequence_7	57.143	14	6	0	42	55	12	25	0.38	13.1
Sequence_1	Sequence_8	61.538	13	5	0	44	56	14	26	0.39	13.1
Sequence_2	Sequence_1	51.852	27	11	1	27	51	26	52	0.030	16.2
Sequence_2	Sequence_1	45.833	24	10	1	35	55	24	0.21	13.9
Sequence_2	Sequence_3	66.667	9	3	0	20	28	19	27	2.9	10.8
Sequence_2	Sequence_4	60.000	5	2	0	44	48	57	61	8.0	 9.6
Sequence_2	Sequence_5	57.143	14	5	1	41	54	56	68	0.68	12.3
Sequence_2	Sequence_7	71.429	7	2	0	11	17	0.65	12.3
Sequence_3	Sequence_1	60.000	5	2	0	18	22	39	43	8.5	 9.2
Sequence_3	Sequence_4	31.250	16	11	0	32	47	46	61	2.9	10.8
Sequence_3	Sequence_5	56.250	16	7	0	16	31	23	38	0.075	15.0
Sequence_3	Sequence_5	50.000	6	3	0	8	13	31	36	3.3	10.4
Sequence_3	Sequence_6	87.500	8	1	0	30	37	61	68	0.061	15.4
Sequence_3	Sequence_7	80.000	5	1	0	58	62	23	27	6.2	 9.6
Sequence_3	Sequence_7	36.364	11	7	0	18	28	16	26	9.1	 9.2
Sequence_3	Sequence_8	41.667	12	7	0	3	14	28	39	0.059	15.4
Sequence_3	Sequence_8	40.000	15	9	0	56	70	17	31	4.5	10.0
Sequence_4	Sequence_2	35.714	14	9	0	57	70	44	57	4.7	10.0
Sequence_4	Sequence_3	35.556	45	23	1	23	61	47	1.0	11.9
Sequence_4	Sequence_5	36.000	25	16	0	13	37	16	40	0.57	12.7
Sequence_4	Sequence_6	100.000	5	0	0	43	47	44	48	0.20	13.9
Sequence_4	Sequence_7	45.455	11	6	0	8	18	12	22	0.58	12.7
Sequence_5	Sequence_10	50.000	8	4	0	58	65	16	1.8	11.2
Sequence_5	Sequence_10	50.000	6	3	0	59	64	30	35	9.3	 9.2
Sequence_5	Sequence_2	57.143	14	5	1	56	68	41	54	0.35	13.1
Sequence_5	Sequence_3	50.000	30	8	1	9	38	31	0.003	18.9
Sequence_5	Sequence_3	50.000	6	3	0	31	36	13	3.7	10.4
Sequence_5	Sequence_4	36.364	33	21	0	4	36	33	0.13	14.2
Sequence_5	Sequence_4	75.000	8	2	0	5	12	61	68	0.18	13.9
Sequence_5	Sequence_6	55.556	9	4	0	3	11	60	68	0.35	13.1
Sequence_5	Sequence_6	42.857	56	26	3	2	57	58	0.44	13.1
Sequence_5	Sequence_9	60.000	10	4	0	42	51	18	2.0	11.2
Sequence_6	Sequence_1	45.455	11	6	0	44	54	19	0.43	13.1
Sequence_6	Sequence_1	80.000	5	1	0	22	26	63	67	5.1	10.0
Sequence_6	Sequence_1	32.558	43	24	1	7	49	15	52	9.2	 9.2
Sequence_6	Sequence_3	87.500	8	1	0	61	68	30	37	0.068	15.0
Sequence_6	Sequence_4	100.000	5	0	0	44	48	43	47	0.27	13.5
Sequence_6	Sequence_4	55.556	9	4	0	25	33	15	3.0	10.8
Sequence_6	Sequence_5	45.238	42	17	3	23	58	16	57	0.097	14.6
Sequence_6	Sequence_7	66.667	6	2	0	10	15	18	23	0.41	13.1
Sequence_6	Sequence_8	55.556	9	4	0	37	45	15	23	9.0	 9.2
Sequence_7	Sequence_10	75.000	8	2	0	32	39	28	35	1.4	11.5
Sequence_7	Sequence_10	66.667	3	1	0	44	46	26	28	8.8	 9.2
Sequence_7	Sequence_2	87.500	8	1	0	54	61	0.11	14.6
Sequence_7	Sequence_2	83.333	6	1	0	2	7	12	17	0.75	12.3
Sequence_7	Sequence_3	85.714	7	1	0	23	29	58	64	0.14	14.2
Sequence_7	Sequence_3	44.444	27	11	1	16	42	18	40	0.33	13.5
Sequence_7	Sequence_3	39.535	43	8	2	32	65	30	63	1.6	11.5
Sequence_7	Sequence_4	35.484	31	13	1	13	43	32	0.50	12.7
Sequence_7	Sequence_5	57.895	19	6	1	28	44	27	45	1.3	11.5
Sequence_7	Sequence_6	66.667	12	2	1	18	27	10	21	0.75	12.3
Sequence_7	Sequence_6	52.941	17	6	1	25	39	52	68	7.1	 9.6
Sequence_8	Sequence_10	60.000	15	4	1	4	18	21	33	0.29	13.5
Sequence_8	Sequence_1	58.333	12	5	0	42	53	34	45	1.6	11.5
Sequence_8	Sequence_1	54.545	11	5	0	14	24	44	54	2.1	11.2
Sequence_8	Sequence_1	50.000	6	3	0	16	21	28	33	8.4	 9.2
Sequence_8	Sequence_2	83.333	6	1	0	48	53	39	44	0.36	13.1
Sequence_8	Sequence_3	44.118	34	17	1	28	59	36	0.008	17.7
Sequence_8	Sequence_3	40.000	15	9	0	17	31	56	70	4.3	10.4
Sequence_8	Sequence_4	71.429	7	2	0	48	54	15	21	0.69	12.3
Sequence_8	Sequence_9	62.500	8	3	0	16	23	62	69	3.3	10.4
Sequence_9	Sequence_10	43.750	16	9	0	53	68	17	32	0.63	12.7
Sequence_9	Sequence_10	60.000	5	2	0	37	41	57	61	9.3	 9.2
Sequence_9	Sequence_1	100.000	3	0	0	41	43	52	54	2.4	10.8
Sequence_9	Sequence_2	38.462	13	8	0	33	45	13	8.0	 9.6
Sequence_9	Sequence_5	41.935	31	15	1	31	58	17	47	0.35	13.1
Sequence_9	Sequence_5	41.667	24	14	0	9	32	42	65	2.2	11.2
Sequence_9	Sequence_5	41.176	17	10	0	33	49	45	61	9.3	 9.2
Sequence_9	Sequence_8	57.143	14	6	0	56	69	10	23	0.43	13.1
Sequence_9	Sequence_8	83.333	6	1	0	31	36	65	70	0.97	11.9
