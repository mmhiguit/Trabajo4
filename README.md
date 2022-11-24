# ANÁLISIS DE LA ABUNDANCIA RELATIVA DE DINUCLEÓTIDOS EN CEPAS PATÓGENAS Y NO PATÓGENAS DEL GÉNERO Pseudomonas.



# CONSIDERACIONES GENERALES

Esta rutina permite realizar la comparación entre la abundancia relativa de dinucleotidos (RDA) a partir de la secuencia del genoma bacteriano,adcionalmente, 
permite determinar si existe o no una diferencia significativa entre el RDA de cepas bacterianas diferentes. En este caso se analizaron secuencias genómicas de
dos cepas doferentes de Pseudomonas, una patógena (Pseudomonas aeruginosa) y otra no patógena (Pseudomonas putida). También se realizo el análisis del RDA en
una sequencia de 6156701 nt, la cual fue generada aleatoriamente con el paquete de python random, la longitud de la secuencia fue elegida, con el fin de que esta
conservara la misma longitud del genoma de las bacterias escogidas. A continuación, se detalla la información de las secuencias empleadas para el análisis.


<img width="694" alt="image" src="https://user-images.githubusercontent.com/116923271/202936152-ad9c6230-aab3-48bb-8414-8313647a17e3.png">


# ABUNDANCIA RELATIVA DE DINUCLEÓTIDOS

La abundancia relativa de dinucleótidos es frecuentemente utilizada para realizar comparaciones entre genomas de diferentes organismos y entre diferentes especies 
bacterianas (Srividhya et al., 2007), esta métrica es altamente conservada entre oraganismos de la misma especie. Por esta razón,  se espera que la RDA permita establecer diferencias entre dos especies de Pseudomonas, una patógena (Pseudomonas aeruginosa) y otra no patógena (Pseudomonas putida). Con el fin de garantizar que esta métrica si permite visualizar estas diferencias, también 
fueron realizados los análisis con una secuencia de la misma longitud generada de manera aleatoria.

Para calcular la abundancia relativa se utiliza la siguiente ecuación:

<img width="504" alt="image" src="https://user-images.githubusercontent.com/116923271/202954394-d6b86c25-c7a8-4d18-8dc0-99e52676af0a.png">


Para calcular los RDA a partir de los archivos fasta de las tres secuencias se utilizó el paquete de Python DinuQ, desarrollado por Lytras et al., 2020
(https://github.com/spyros-lytras/dinuq/blob/master/README.md). Este programa permite generar un diccionario, con las abundancias relativas de cada dinucleotido
para cada una de las secuencias analizadas. 







# RUTINA PARA EL ANÁLISIS DE LA ABUNDANCIA RELATIVA DE LOS DINUCLEÓTIDOS

Con esta rutina se espera determinar si existen diferencias significativas entre la abundancia relativa de los dinucleótidos de la cepa patógena y la 
no patógena.


# Paquetes requeridos

Para el funcionamiento de esta rutina será necesario la instalación de los siguientes paquetes.

		import os
		import sys
		import pandas as pd
		import matplotlib.pyplot as plt
		import seaborn as sns
		import dinuq as dn # This package provides a range of metrics for quantifying dinucleotide in genetic sequences. (pip install dinuq)
		from scipy.stats import f_oneway
		import random
		import scipy.stats as stats
		from scipy.stats import tukey_hsd
		from statsmodels.stats.multicomp import pairwise_tukeyhsd


# Generación de la secuencia aleatoria

Con el fin de determinar si esta métrica permitía detectar diferencias entre secuencias diferentes y los dinucleótidos que generaban dichas diferencias, se creó una secuencia aleatoria con el paquete random. La función de esta secuencia es actuar como control negativo en los análisis que se realizarán, permitiendo observar si existen o no patrones conservadps en los perfiles generados por los valores de RDA para las secuencias de los genomas bacterianos.


# Cálculo de los RDA

Como se mencionó previamente, para calcular los valores RDA de cada uno de los 16 dinucleótidos se utilizó el programa dinuq y se ejecutaron los siguientes
comandos.




					#Create random sequence with a length of 6156701 to compare with pathogen and No pathogen sequences

					 print(f'Create random sequence')
					 length= 6156701 #Define length of sequence

					 bases=["A","G","C","T"] #supply list of nucleotides

					 #Generate the random sequence

				   from random import choice
					#create empty sequence
					 random_sequence = ""


					#randomly select base and add to sequence(repeat 100 times for length=6156701)
						for i in range(length):
						base=choice(bases)
						random_sequence+=base

Después de tener la secuencia, es necesario que esta se guarde en un archivo fasta, el cual, será requerido más adelante por el cáculo de los RDA
con el programa DinuQ. 

					#define function to create a fasta file with the random sequence

					f=open("random_seq.fasta","w")
					f.write(f'>random \n{random_sequence}')

					# print(sequence) # Just to verify sequence

Este archivo se creará con el nombre 'random_seq.fasta' en la carpeta donde se encuentre el código de la rutina.

#Cálculo de los RDAs paras las tres secuencias (Patógena, No patógena y Random)

Para calcular los RDAs desde la consola de python, es necesario definir las tres secuencias, de la siguiente manera.


      #Activate to run code python console
			# Pathogen = 'Pseudomonas_aeruginosa.fasta'
			# No_pathogen= 'Pseudomonas_putida.fasta'
			# Other ='random_seq.fasta'
			
Por otro lado, si se quiere ejecutar el código desde la terminal, será necesario inactivar la instrucción anterior y correr el comando sys.argv[], 
de la siguiente manera:

			#Activate to run code from terminal and use: python iRutina_P4.py Pseudomonas_aeruginosa.fasta Pseudomonas_putida.fasta random_seq.fasta
			Pathogen = sys.argv[1]
			No_pathogen = sys.argv[2]
			Other = sys.argv[3]
			
Si se utiliza esta última opción, el código podrá correrse desde la terminal empleando el siguiente input: 
 <% python Pseudomonas_aeruginosa.fasta Pseudomonas_putida.fasta random_seq.fasta>

<img width="1352" alt="image" src="https://user-images.githubusercontent.com/116923271/202961567-45570f81-797a-447b-81a7-d1cfdfd2dc72.png">

# RDAs

Para poder calcular los RDAs, es necesario generar una lista de las 16 combinaciones de dinucleótidos posibles, como se muestra a continuación,

             dinucl_list = ['TpA','TpC','TpG','TpT','ApA',
               'ApC','ApG','ApT','CpA','CpC',
              'CpG','CpT','GpA','GpC','GpG','GpT'] #list of dinucleotides
							
Después de generar la lista, se procede a calcular los valores RDA para cada secuencia y sus respectivos data frames


# RDA Cepa Patógena (Pseudomonas aeruginosa)

             
						  rda_pathogen = dn.RDA(Pathogen, dinucl_list)
							print(rda_pathogen)


							RDA1_df = pd.DataFrame.from_dict(rda_pathogen).reset_index() # Create data frame with RDAs values
							RDA1_df.rename(columns = {'index':'Dinucleotide'}, inplace = True)
							print(RDA1_df)
<img width="1383" alt="image" src="https://user-images.githubusercontent.com/116923271/202963011-97986e9d-5511-4f15-a666-9e86d2821d8b.png">


# RDA cepa No patógena (Pseudomonas_putida.fasta)

							rda_nopathogen = dn.RDA(No_pathogen, dinucl_list)
							print(rda_nopathogen)

							RDA2_df = pd.DataFrame.from_dict(rda_nopathogen).reset_index() # Create data frame with RDAs values
							RDA2_df.rename(columns = {'index':'Dinucleotide'}, inplace = True)
							print(RDA2_df)

<img width="1383" alt="image" src="https://user-images.githubusercontent.com/116923271/202963037-64bef1ef-c973-4871-a60d-e3c16483e0fd.png">


# RDA Secuencia control (random_seq.fasta)

							rda_other = dn.RDA(Other, dinucl_list)
							print(rda_other)

							RDA3_df = pd.DataFrame.from_dict(rda_other).reset_index() # Create data frame with RDAs values
							RDA3_df.rename(columns = {'index':'Dinucleotide'}, inplace = True)
							print(RDA3_df)
							
<img width="1383" alt="image" src="https://user-images.githubusercontent.com/116923271/202963103-047786f8-2f3d-4351-997b-8463433bcf6e.png">


Como se puede observar, ya se generaron los data frame con los resultados del los RDAs de las secuencias analizadas, sin embargo, estos se encuentran en 
diferentes data frame, y es necesario concatenarlos en uno solo. Para unir los tres data frame, primero se unen el df del RDA1 y el RDA 2, luego a este nuevo 
df concatenado, se une el RDA3, como se indica a continuación:

# Unión de los df de RDA1 y RDA2

				# Join RDA's Tables
				df_conc = RDA1_df.join(RDA2_df[list(rda_nopathogen.keys())[0]]) # Concatenate RDA1 and RDA2 data frame
				print(df_conc)

<img width="372" alt="image" src="https://user-images.githubusercontent.com/116923271/202963583-2ad340d3-6d0a-448a-a53a-8a799a7d861e.png">

# Unión del df concatenado de RDA1 y RDA2 (df_conc) con el df de RDA3

				df_conc_2 = df_conc.join(RDA3_df[list(rda_other.keys())[0]]) # Concatenate df_conc and RDA3 data frame
				print(df_conc_2)

<img width="426" alt="image" src="https://user-images.githubusercontent.com/116923271/202964070-460d1a12-c10f-4c08-b66f-2cd6f856d64e.png">


Finalmente, se realiza la edición de la tabla para crear las categorías y poder realizar la elaboración de los gráficos y los análisis estadísticos

					print(f'editing concatening data frame')
					# Add names to columns in the fd_conc_2
					dfm = df_conc_2.melt('Dinucleotide', var_name='Accession', value_name='RDA')
					print(dfm)
					
<img width="321" alt="image" src="https://user-images.githubusercontent.com/116923271/202964471-89d2891c-90e6-4f07-b024-dbc648736a82.png">


# Elaboración de gráficos

Con los resultados obtenidos se procede a elaborar seis gráficos diferentes, barplot, relplot, scatteplot, boxplot, y lineplot,  utilizando el paquete seaborn de python (Fig.1). De los gráficos generados, el que mejor ilustra los valores extremos en los valores RDA correspondiente a los tres organismos, es el lineplot (Fig.2). En este se evidencian claramente las diferencias entre los perfiles de RDA de la secuencia patógena y la secuencia no patógena. Estas diferenncias son generadas principalmente por la abundancia relativa de los dinucleótidos AT, CA, CG, TG y TC, como se observa en la fig 2. 
que esta métrica resulta es adecuada para evidenciar cambios en los genomas bacterianos.

  # Figura 1: 
  Gráficos para ilustrar la abundancia relativa en los dinucleótidos de bacterias patógenas (Pseudomonas aeruginosa) y no patógenas (Pseudomonas putida).
  
<img width="978" alt="image" src="https://user-images.githubusercontent.com/116923271/203669386-8aebb2cd-18bc-4041-ac09-cfa878672a65.png">


  # Figura 2:  
  Lineplot de la abundancia relativa en los dinucleótidos de bacterias patógenas (Pseudomonas aeruginosa) y no patógenas (Pseudomonas putida).
  
<img width="559" alt="image" src="https://user-images.githubusercontent.com/116923271/203671365-48d3da74-374f-4446-acb9-b5506a1d2607.png">




# Análisis estadístico

El análisis estadístico se realizó considerando un alfa = 0.05 y bajo el siguiente juego de hipotesis.

Ho:Las medias entre los RDAs son iguales, no hay diferencia significativa
Ha:Las medias entre los RDAs son diferentes, hay diferencia significativa

Dada estas hipotesis se evaluará lo siguiente, si el valor p es mayor a 0.05 no se cuenta con suficiente evidencia para rechazar Ho y por ende se concluye 
que no hay diferencias significativas, por otro lado, si el valor p es menor o igual a 0.05 es posible rechazar Ho y puede concluirse que hay diferencias 
significativas.

# ANÁLISIS DE LA VARIANZA - ANOVA
El anova se realizó utilizando la herramieta stats del paquete scipy.stats, empleando el siguiente comando.

		fvalue, pvalue = stats.f_oneway(df_conc_2['NC_002516.2'], df_conc_2['NC_021505.1'], df_conc_2['random'])
		print(fvalue, pvalue)

Obteniendo el siguiente resultado

<img width="342" alt="image" src="https://user-images.githubusercontent.com/116923271/202966424-dfea50e0-b49d-4442-915b-c431403d7453.png">

Dado que el valor p fue igual a 0.9999 y es mayor a 0.05, se concluye que no se cuenta con evidencia para rechazar la hipotesis nula, indicando que las medias son iguales, y por tanto no hay diferencias significativas entre los RDAs de los organismos analizados. Sin embargo, en los gráficos se evidenció que hay diferencias entre los perfiles de RDA de las secuencias, esta incongruencia, se debe principalmente a que el ANOVA no tiene el poder para detectar 
las diferencias presentes en este tipo de análisis, en los cuales, es recomendado realizar el análisis estadístico Least significant Difference (LSD) (Wang et al.,2022), sin embargo, al realizar busqueda en diferentes páginas de programación y en espacios de discusión de python, no se encontraron paquetes que permitieran realizarlo en python y hasta el momento solo puede calcularse usando R studio. Por esta razón y con el objetivo de poder observar las diferencias entre los pares de secuencias analizadas se llevo a cabo el test Tukey HSD (Honestly-significant-difference) disponible en el paquete scipy.stats. Este test permite hacer comparaciones multiples al analizar las diferencias en pareja entre las medias, permitiendo una mayor resolución en el análisis de los valores de RDA entre pares de secuencias. Para realizar este análisis se corrió en puthon lo siguiente,

# Test HSD (Honestly-significant-difference) de Tukey.
			print(f'Tukey_HSD')
			Tukey_HSD = tukey_hsd(df_conc_2['NC_002516.2'], df_conc_2['NC_021505.1'], df_conc_2['random'])

			print(Tukey_HSD)

Obteniendo el siguiente resultado,

  <img width="537" alt="image" src="https://user-images.githubusercontent.com/116923271/203675645-2c929bdf-cc1b-4450-9656-ed27340dd6ef.png">


Para analizar el anterior resultado se debe tener presente que:

0 = NC_002516.2 (Pathogen)
1 = NC_021505.1 (No pathogen)
2 = random  (Random_seq)

Los datos obtenidos en el Tukeys'HSD, confirman que no existen diferencias significativas entre las medias de los RDAs entre las bacterias analizadas, por tal motivo, y considerando que estos resultados no coinciden con lo observado en las gráficas se sugiere realizar el least significant differences, que permita detectar las diferencias, entre el organismo patógeno y no patógeno, en la abundadncia relativa de dinucleótidos AT, CA, CG, CT, TG y TC. También, se sugiere realizar este análisis con un mayor número de genomas bacterianos y comparando las diferencias entre cada uno de los dinucleótidos que presentaron valores extremos en este estudio.




# BIBLIOGRAFÍA


Srividhya, K. V., Alaguraj, V., Poornima, G., Kumar, D., Singh, G. P., Raghavenderan, L., ... & Krishnaswamy, S. (2007). Identification of prophages in bacterial genomes by dinucleotide 
relative abundance difference. PLoS One, 2(11), e1193.

Wang, Q., Lyu, X., Cheng, J., Fu, Y., Lin, Y., Abdoulaye, A. H., ... & Xie, J. (2022). Codon usage provides insights into the adaptive evolution of mycoviruses in their associated 
fungi host.International journal of molecular sciences, 23(13), 7441.


