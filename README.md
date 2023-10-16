# Estimation of the localization space of amino acid residues 101 and 195 of histone [H1.0](https://www.rcsb.org/structure/6hq1) 

### Task: Evaluating FRET Efficiency Values for Nucleosome with Linkers
In this task, we use the experimental FRET efficiency values for a system consisting of a nucleosome with a linker sequence. Where the fluorophore Cy3 is introduced into the DNA sequence, in the linker region of the nucleosome at positions 10-15-20-25-25-30-35 nucleotide park (n.p.). Spacers from two different companies, Syntol© and Lumiprobe©, were used to collect experimental data for the Cy3 fluorophore. 

 - [Syntol©](https://github.com/NVKristovs/nucl_FRET_analysis/blob/main/images/Syntol_Lumiprobe_Cy3.jpg) linkers at positions 15, 20 had a spacer length of 16 Å. 
 - [Lumiprobe©](https://github.com/NVKristovs/nucl_FRET_analysis/blob/main/images/Syntol_Lumiprobe_Cy3.jpg) linkers at positions 10, 25, 30, 35 had a spacer length of 24 Å. 
<img src=" https://github.com/NVKristovs/nucl_FRET_analysis/blob/main/images/Syntol_Lumiprobe_Cy3.jpg" width="600" height="800" style="max-width: 100%;">
The fluorophore [Су5](https://github.com/NVKristovs/nucl_FRET_analysis) was introduced into the histone H1 sequence at positions 101 and 195 amino acids, with a spacer length of 11 Å.
Reaction with [Maleimide](https://ru.lumiprobe.com/protocols/protein-maleimide-labeling) was used to label the protein.

The library [LabelLib](https://github.com/Fluorescence-Tools/LabelLib) was used to determine the available amount of Cy3 and the limitations of the localization region of Cy5 introduced into the protein sequence of 101 and 195 n.b. hrs. This library creates a three-dimensional matrix describing the label localization region based on the specified parameters.

Using the LabelLib library, we determined the regions of available volume of Cy3 and possible localization of Cy5. By analyzing the coordinates of the available volume voxels obtained from LabelLib for Cy3, we determined which voxels in the Cy5 localization region allow FRET efficiency values that fall within the range of one standard deviation of FRET efficiency from experimental data. Thus, we limited the volume in which the Cy5 fluorophore amino acid could be located. We retained only those voxel coordinates that describe the experimental FRET efficiency values. We also calculated the relative probability of finding the Cy5 fluorophore in a particular voxel. This probability is calculated as the sum of the number of voxel hits within the range of one standard deviation of FRET efficiency for a fluorophore introduced at one of the positions (10/15/20/25/30/35) divided by the maximum number of hits for that Cy3 position and summed for all 6 Cy3 introduction positions.

#### Sources:

- Dimura, M., Peulen, T.O., Hanke, C.A., Prakash, A., Gohlke, H. and Seidel, C.A., 2016. [Quantitative FRET studies and integrative modeling unravel the structure and dynamics of biomolecular systems. Current opinion in structural biology](https://www.sciencedirect.com/science/article/pii/S0959440X1630197X?via%3Dihub), 40, pp.163-185.
--------------------------------------------------------------------------------
# Оценка пространства локализации остатков аминокислот 101 и 195 гистона [Н1.0](https://www.rcsb.org/structure/6hq1) 

### Задача: Оценить значения эффективности FRET для системы: нуклеосома с линкерами. 
В данной задаче мы используем эксперементальные значения эффективность FRET для системы, состоящей из нуклеосомы с линкерной последовательностью. Где флуорофор Су3 введен в последовательность ДНК, в линкерной области нуклеосомы в позиции 10-15-20-25-30-35 нуклеотидных парк (н.п.). При сборе эксперементальных данных для флуорофора Cy3 использовались спейсеры от двух различных компаний: Syntol© и Lumiprobe©. 
Где линкеры от компании [Syntol©](https://github.com/NVKristovs/nucl_FRET_analysis/blob/main/images/Syntol_Lumiprobe_Cy3.jpg) в позициях 15, 20. — длина спейсера 16 Å.
Линкеры от компании [Lumiprobe©](https://github.com/NVKristovs/nucl_FRET_analysis/blob/main/images/Syntol_Lumiprobe_Cy3.jpg) в позициях 10, 25, 30, 35. — длина спейсера 24 Å.

Флуорофор [Су5](https://github.com/NVKristovs/nucl_FRET_analysis) введен в последовательность гистона Н1 в позиции 101 и 195 аминокислот (а.к.), длина спейсера 11 Å.

Для маркировки белка использовалась реакция с [Малеимидом](https://ru.lumiprobe.com/protocols/protein-maleimide-labeling)

Для определения доступного объема Cy3 и ограничений области локализации Сy5 введенного в белковую последовательность чреез 101 и 195 н.к. была использована библиотека [LabelLib](https://github.com/Fluorescence-Tools/LabelLib). Эта библиотека создает трехмерную матрицу, описывающую область локализации метки на основе заданных параметров.

С использованием библиотеки LabelLib мы определили области доступного объема Cy3 и возможной локализации Cy5. Путем анализа координат вокселей доступного объема, полученных из LabelLib для Cy3, мы определили, какие воксели в области локализации Cy5 позволяют получить значения эффективности FRET, входящие в диапазон одного стандартного отклонения эффективности ФРПЭ из экспериментальных данных. Таким образом, мы ограничили объем, в котором может находиться аминокислота с флуорофором Cy5. Мы сохранили только те координаты вокселей, которые описывают экспериментальные значения эффективности FRET. Также была рассчитана относительная вероятность нахождения флуорофора Cy5 в конкретном вокселе. Эта вероятность вычисляется как сумма количества попаданий вокселя в диапазон одного стандартного отклонения эффективности ФРПЭ для флуорофора, введенного в одну из позиций (10/15/20/25/30/35), деленная на максимальное количество попаданий для данной позиции  Cy3 и суммированная для всех 6 позиций введения Cy3.

#### Источники:

- Dimura, M., Peulen, T.O., Hanke, C.A., Prakash, A., Gohlke, H. and Seidel, C.A., 2016. [Quantitative FRET studies and integrative modeling unravel the structure and dynamics of biomolecular systems. Current opinion in structural biology](https://www.sciencedirect.com/science/article/pii/S0959440X1630197X?via%3Dihub), 40, pp.163-185.


### Результаты 
1.  
2. 
3. Полученные пространства являются грубой оценкой локализации искомых аминокислот и согласуются с современными представлениями о структуре Н1 гистона.
