//  Segundo Trabalho de MNC 
//  Sistemas Lineares
// Leonardo Oliveira

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <locale.h>

#define max 10

double Determinante (int ordem, double matriz[][max]); 

int Menu(){
	
	int x; 
	
	printf("\nEscolha a opção: "); 
	printf("\n\n1 - Cálculo de Determinante"); 
	printf("\n10 - Fechar programa"); 
	printf("\n\nOpção: "); 
	fflush(stdin); 
	scanf("%d",&x); 
	
	return x; 
}

void MenorPrincipal (double matriz[max][max],double Menor[max][max], int linha, int coluna, int ordem){
	
	int aux_l = 0, aux_c = 0; 
	
	for (int i=0; i<ordem; i++){
		if(i == linha) continue; 
		for (int j=0; j<ordem; j++){
			if(j == coluna)	continue; 
			Menor[aux_l][aux_c++] = matriz[i][j]; 			
		}
		aux_c = 0; 
		aux_l++; 
	}	
}

void MelhorCaminho (double matriz[max][max],int ordem, int *lin_col, int *n){
	
	int cont_i, cont_j, melh_i, melh_j, lin[ordem], col[ordem]; 
	cont_i = cont_j = melh_i = melh_j = (*lin_col) = 0;	 
	
	for (int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++)
			if(matriz[i][j] == 0) cont_i++; 
		lin[i] = cont_i; 
		cont_i = 0; }
	
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++)
			if(matriz[j][i] == 0) cont_j++;
		col[i] = cont_j; 
		cont_j = 0; }
	
	for(int i=0; i<ordem; i++){
		if(lin[i] != 0)
			if(lin[melh_i] < lin[i])
				melh_i = i; }
	
	for(int i=0; i<ordem; i++){
		if(col[i] != 0)
			if(col[melh_j] < col[i])
				melh_j = i; }
	
	if(lin[melh_i] >= col[melh_j]){
		(*lin_col) = 0; 
		(*n) = melh_i; 
	}else{
		(*lin_col) = 1; 
		(*n) = melh_j; 
	}	
}

void MostraMatriz ( int ordem, double matriz[][max]){	
	for(int i=0; i<ordem; i++){
		printf("\n"); 
		for(int j=0; j<ordem; j++)
			printf(" %lf ",matriz[i][j]); 
	}
}

double Cofator (int ordem, int i, int j, double matriz[][max]){
	return pow((-1),i+j) * Determinante(ordem-1,matriz);  
}

double Determinante (int ordem, double matriz[][max]){
	
	int i, j, lin_col, num;
	double maux [max][max]; 
	int soma = 0; 
	
	if(ordem == 2) return (matriz[0][0] * matriz[1][1] - matriz[0][1]*matriz[1][0]);		
	
	if(ordem == 1) return matriz[0][0];
	
	MelhorCaminho(matriz,ordem,&lin_col,&num);
	
	switch(lin_col){
		case 0: 	i = num; 
					for(j=0; j<ordem; j++){
					
						MenorPrincipal (matriz,maux,i,j,ordem); 
						soma += matriz[i][j] * Cofator (ordem,i,j,maux);
					}
					break; 
		case 1: 	j = num; 
					for(i=0; i<ordem; i++){
					
						MenorPrincipal (matriz,maux,i,j,ordem); 
						soma += matriz[i][j] * Cofator (ordem,i,j,maux);
					}
					break; 	}
	
	return soma; 
}

LeMatriz(double matriz[][max], int *ordem){
	
	printf("\nDigite a ordem da matriz: "); 
	scanf("%d",&(*ordem)); 
	printf("\nDigite a matriz: "); 
	for(int i=0; i<(*ordem); i++)
		for (int j=0; j<(*ordem); j++)
			scanf("%lf",&matriz[i][j]); 
	
}

int main(){

	setlocale(LC_ALL,"Portuguese");
	int ordem, det, aux = 0, num, op; 
	double matriz[max][max], menor[max][max]; 
	
	//freopen("matriz.txt","r",stdin); 
	
	do{
		op = Menu(); 
		switch(op){
			case 1: 	LeMatriz(matriz,&ordem); 
						printf("\nMatriz: \n"); 
						MostraMatriz(ordem, matriz); 
						printf("\n\nDet(A) = %d",Determinante(ordem,matriz));
						fflush(stdin); 
						system("pause"); 
						break; 
			case 2: 	
						break; 
		}
	}while(true); 
	
	scanf("%d",&ordem); 
	for(int i=0; i<ordem; i++)
		for(int j=0; j<ordem; j++)
			scanf("%lf",&matriz[i][j]); 
			
	
	printf("\nOrdem: %d",ordem); 
	MostraMatriz(ordem,matriz); 
	
	/*printf("\nMelhor caminho para calcular a determinante: "); 
	MelhorCaminho(matriz,ordem,&aux,&num); 
	
	if(!aux)
		printf("\nLinha de numero %d",num+1); 
	else printf("\nColuna de numero %d",num+1); 
	
	printf("\n\n");*/
	
	det = Determinante(ordem,matriz); 
	
	printf("\nDeterminante: %d",det); 
	
	printf("\nFim do programa"); 
}

