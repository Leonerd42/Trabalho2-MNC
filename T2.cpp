/*
Alunos: 
Lucas Henrique Russo do Nascimento, RA: 171025709
Bruna Lika Tamake, RA: 171024427
Barbara Este Fernandez, RA:161025901
Leonardo Silva de Oliveira, RA: 171025903
*/
#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include<ctype.h>
#include<math.h>
#define MAX 10
int n;
//------------------------------------------
//Determinante:
float **ReduzMatriz(int ordem, int i, float **matriz, float **matrizReduzida){
	/*int j, k;
	float **p=matrizReduzida, *p2;
	for(j = 1; j < ordem; j++){
		p2 = *p++;
		for(k = 0; k < ordem; k++){
			if(j==i){
				continue;
			}
			else{
				*p2++ = matriz[j][k];
			}
		}
	}
	return matrizReduzida;*/
}
float Determinante(int ordem, float **matriz){
	/*float det=0;
	if(ordem==1){
		det=matriz[1][1];
		return det;
	}
	int cofator;
	float **matrizReduzida;
	for(int i = 0; i < ordem; i++){
		if(i%2==1){
			ReduzMatriz(ordem, i, matriz, matrizReduzida);
			cofator=(-1)*Determinante(ordem-1, matrizReduzida);
		}
		else{
			ReduzMatriz(ordem, i, matriz, matrizReduzida);
			cofator=1*Determinante(ordem-1, matrizReduzida);
		}
		det += matriz[0][i]*cofator;
	}*/
}
//------------------------------------------
int menu_metodos(){
	int op;
	do{
		system("cls");
		printf("\n\tMENU\t\n");
		printf("\n\t1-Determinante \n\t2-Sistema Triangular Inferior \n\t3-Sistema Triangular Superior \n\t4-Decomposicao LU \n\t5-Cholesky \n\t6-Gauss Compacto \n\t7-Gauss Jordan \n\t8-GPP SEM troca de linhas \n\t9-Jacobi \n\t10-Gauss Seidel \n\t11-Matriz Inversa\n\t12-Fechar programa\n");
		printf("\n\tDigite a opcao:");
		scanf("%d",&op);
		if(op<1 || op>12){
			printf("\n\nOpcao Invalida!\n");
		}
	}while(op<1 || op>12);
	return op;
}
//------------------------------------------
void ColetaDados(int *ordem, int opcao){
	if(opcao==1){ // coleta dados para calculo do determinante
		system("cls");
		printf("\n\tDigite a ordem da matriz:");
		scanf("%d", &*ordem);
	}
}
//------------------------------------------
void LerMatriz(float **matriz, int ordem){
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++){
			printf("a%d%d = ",i+1, j+1);
			scanf("%f", &matriz[i][j]);
		}
	}
}
void MostraMatriz(float **matriz, int ordem){
	printf("\n\tMatriz:\n");
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++){
			printf("\t%.2f", matriz[i][j]);
		}
		printf("\n");
	}
}
float **alocarMatriz(int ordem){
	float **m;
	m = (float**)calloc(ordem, sizeof(float*));
	if(m == NULL){
		exit(1);
	}
	for(int i=0; i<ordem; i++){
		*(m+i) = (float*)calloc(n, sizeof(float));
		if(*(m+i) == NULL){
			exit(1);
		}
	}
	return m;
}
//------------------------------------------
int main(){
	int opcao, ordem;
	float det, **matriz; 
	opcao=menu_metodos();
	switch(opcao){
		case 1:
			ColetaDados(&ordem, 1);
			matriz= alocarMatriz(ordem);
			LerMatriz(matriz, ordem);
			//MostraMatriz(matriz, ordem);
			det=Determinante(ordem, matriz);
			printf("det(A) = %f", det);
			break;
		/*case 2:
			break;
		case 3:
			break;
		case 4:
			break;
		case 5:
			break;
		case 6:
			break;
		case 7:
			break;
		case 8:
			break;
		case 9:
			break;
		default:
			return 0;
			break;*/
	}
}
