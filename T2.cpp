/***************************************************
		Trabalho de MNC - Sistemas Lineares
Lucas Henrique Russo do Nascimento, RA: 171025709
Bruna Lika Tamake, RA: 171024427
Barbara Este Fernandez, RA:161025901
Leonardo Silva de Oliveira, RA: 171025903

*****************************************************/
#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include<ctype.h>
#include<math.h>
#define max 10

double Determinante (int ordem, double matriz[][max]);

//  CALCULO DE DETERMINANTE

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

//  FIM DAS ROTINAS DO CALCULO DE DETERMINANTE

/*************************************************************
		     Cálculo dos Sistemas Triangulares
**************************************************************/

void SistemaTriangularSuperior(int ordem, double matriz[][max], double b[], double x[]){ //okay, funciona!!
	
	double somatorio;  
	
	for(int i=0; i<ordem; i++){
		somatorio = 0; 
		for (int j=0; j<i; j++){
			somatorio += x[ordem-1-j] * matriz[ordem-1-i][ordem-1-j]; 
		}
		x[ordem-1-i] = (b[ordem-1-i] - somatorio) / matriz[ordem-1-i][ordem-1-i]; 
	}
}

void SistemaTriangularInferior(int ordem, double matriz[][max], double b[], double x[]){
	
	double somatorio; 
	
	for(int i=0; i<ordem; i++){
		somatorio = 0; 
		for(int j=0; j<i; j++){
			somatorio += (x[j] * matriz[i][ordem-1-j]); 
		}
		x[i] = (b[i] - somatorio) / matriz[i][ordem-1-i];
	}
}



int MenuMetodos(){
	int op = 0;
	do{
		system("cls");
		fflush(stdin);
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
void ColetaDados(int *ordem, int opcao, double matriz[][max]){
	if(opcao==1){ // coleta dados para calculo do determinante
		system("cls");
		printf("\n\tDigite a ordem da matriz:");
		scanf("%d", &(*ordem));
		printf("\nDigite a Matriz: \n"); 
		for(int i=0; i<(*ordem); i++)
			for(int j=0; j<(*ordem); j++)
				scanf("%lf",&matriz[i][j]);
	}
}

void MostraMatriz(double matriz[max][max], int ordem){
	printf("\n\tMatriz:\n");
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++){
			printf("\t%.2lf", matriz[i][j]);
		}
		printf("\n");
	}
}
//------------------------------------------
int main(){
	int opcao, ordem;
	double det, matriz[max][max]; 	
	double b[max], solucao[max]; 
	
	freopen("matriz.txt","r",stdin); 	//Matriz inferior
	//freopen("matriz2.txt","r",stdin); // Matriz superior
	scanf("%d",&ordem); 
	for(int i=0; i<ordem; i++)
		for(int j=0; j<ordem; j++)
			scanf("%lf",&matriz[i][j]); 
			
	for(int i=0; i<ordem; i++)
		scanf("%lf",&b[i]); 
	
	MostraMatriz(matriz,ordem); 
	
	printf("\nTermos independentes: "); 
	for(int i=0; i<ordem; i++)
		printf("  %lf ",b[i]);
		
	SistemaTriangularInferior(ordem,matriz,b,solucao); 
	//SistemaTriangularSuperior(ordem,matriz,b,solucao);
	
	printf("\nSolução: "); 
	for(int i=0; i<ordem; i++)
		printf("  %lf ",solucao[i]);	
	
	//Comentei o menu pra poder testar as rotinas individualmente 
	//pra facilitar
	/*do{
		opcao = MenuMetodos();
		switch(opcao){
			case 1:		ColetaDados(&ordem, 1, matriz);
						MostraMatriz(matriz,ordem);
						det=Determinante(ordem, matriz); 
						printf("det(A) = %lf\n\t", det);
						system("pause"); 
						break;
			case 2:    
						break; 
		}
	}while(opcao != 12);*/
}
