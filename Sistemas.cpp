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
	return pow((-1),i+j) * 3;  
}

double Determinante (int ordem, double matriz[][max]){
	
	int i, j; 
	
	//for(i=0; i<ordem; i++)
	double maux [max][max]; 
	
	
	
	
}

int main(){

	setlocale(LC_ALL,"Portuguese");
	int ordem; 
	double matriz[max][max];	
	
	freopen("matriz.txt","r",stdin); 
	
	scanf("%d",&ordem); 
	for(int i=0; i<ordem; i++)
		for(int j=0; j<ordem; j++)
			scanf("%lf",&matriz[i][j]); 
			
	
	printf("\nOrdem: %d",ordem); 
	MostraMatriz(ordem,matriz); 
	
	printf("\n\nValor teste: %lf",Cofator(1,2,2,matriz)); 
	
	printf("\nFim do programa"); 
}

