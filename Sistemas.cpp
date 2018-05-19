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


int menu(){
	
	int x; 
	
	printf("\nEscolha a opção: "); 
	printf("\n\n1 - blabla "); 
	printf("\n10 - Fechar programa"); 
	printf("\n\nOpção: "); 
	scanf("%d",&x); 
	
	return x; 
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
	
	int i, j; 
	double maux [max][max]; 
	double soma = 0; 
	
	if(ordem == 2){
		return (matriz[0][0] * matriz[1][1] - matriz[0][1]*matriz[1][0]);		
	}
	
	if(ordem == 1){
		return matriz[0][0];
	}
	
	//linha se mantem fixa 
	i = 0; 
		
	for(j=0; j<ordem; j++){
		if(j == ordem-1){
			
			
		}else{
		
			MenorPrincipal(matriz,maux,i,j,ordem); 
			
		}
		
		//return 
		
		
		//Calcula menorprincipal de Aij
		// um valor * ( -1^i+j * menor principal) 
	}
		
	
}

void MenorPrincipal (double matriz[max][max],double Menor[max][max], int linha, int coluna, int ordem){
	
	int aux_l = 0; 
	int aux_c = 0; 
	
	for (int i=0; i<ordem; i++){
		if(i == linha)
			continue; 
		for (int j=0; j<ordem; j++){
			if(j == coluna)
				continue; 
			Menor[aux_l][aux_c++] = matriz[i][j]; 			
		}
		aux_c = 0; 
		aux_l++; 
	}	
}

int main(){

	setlocale(LC_ALL,"Portuguese");
	int ordem; 
	double matriz[max][max];	
	double menor[max][max]; 
	
	freopen("matriz.txt","r",stdin); 
	
	scanf("%d",&ordem); 
	for(int i=0; i<ordem; i++)
		for(int j=0; j<ordem; j++)
			scanf("%lf",&matriz[i][j]); 
			
	
	printf("\nOrdem: %d",ordem); 
	MostraMatriz(ordem,matriz); 
	
	printf("\nMenor Principal: \n"); 
	MenorPrincipal(matriz,menor,0,2,ordem);
	MostraMatriz(ordem-1,menor);
	
	//printf("\n\nValor teste: %lf",Cofator(1,2,2,matriz)); 
	
	printf("\nFim do programa"); 
}

