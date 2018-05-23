/***************************************************
		Trabalho de MNC - Sistemas Lineares
Lucas Henrique Russo do Nascimento, RA: 171025709
Bruna Lika Tamake, RA: 171024427
Leonardo Silva de Oliveira, RA: 171025903

*****************************************************/
#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include<ctype.h>
#include<math.h>
#define max 10

void MostraMatriz(double matriz[max][max], int ordem);
double Determinante (int ordem, double matriz[][max]);

/*************************************************************
		     Cálculo dos Determinantes
**************************************************************/
void MenorPrincipal (double matriz[max][max],double Menor[max][max], int linha, int coluna, int ordem){
	int aux_l = 0, aux_c = 0; 
	for (int i=0; i<ordem; i++){
		if(i == linha){
			continue; 
		}
		for (int j=0; j<ordem; j++){
			if(j == coluna)	{
				continue; 
			}
			Menor[aux_l][aux_c++] = matriz[i][j]; 			
		}
		aux_c = 0; 
		aux_l++; 
	}	
}
//-----------------------------------------------------------
void MelhorCaminho (double matriz[max][max],int ordem, int *lin_col, int *n){
	int cont_i, cont_j, melh_i, melh_j, lin[ordem], col[ordem]; 
	cont_i = cont_j = melh_i = melh_j = (*lin_col) = 0;	 
	for (int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++){
			if(matriz[i][j] == 0) {
				cont_i++; 
			}
		}
		lin[i] = cont_i; 
		cont_i = 0; 
	}
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++){
			if(matriz[j][i] == 0){
				cont_j++;
			}
		}
		col[i] = cont_j; 
		cont_j = 0; 
	}
	for(int i=0; i<ordem; i++){
		if(lin[i] != 0){
			if(lin[melh_i] < lin[i]){
				melh_i = i; 
			}
		}
	}
	
	for(int i=0; i<ordem; i++){
		if(col[i] != 0){
			if(col[melh_j] < col[i]){
				melh_j = i; 
			}
		}
	}
	if(lin[melh_i] >= col[melh_j]){
		(*lin_col) = 0; 
		(*n) = melh_i; 
	}
	else{
		(*lin_col) = 1; 
		(*n) = melh_j; 
	}	
}
//-----------------------------------------------------------
double Cofator (int ordem, int i, int j, double matriz[][max]){
	return pow((-1),i+j) * Determinante(ordem-1,matriz);  
}
//-----------------------------------------------------------
double Determinante (int ordem, double matriz[][max]){
	int i, j, lin_col, num, soma=0;
	double maux [max][max]; 
	if(ordem == 2){ 
		return (matriz[0][0] * matriz[1][1] - matriz[0][1]*matriz[1][0]);		
	}
	if(ordem == 1) {
		return matriz[0][0];
	}
	MelhorCaminho(matriz,ordem,&lin_col,&num);
	switch(lin_col){
		case 0: 	
			i = num; 
			for(j=0; j<ordem; j++){
				MenorPrincipal (matriz,maux,i,j,ordem); 
				soma += matriz[i][j] * Cofator (ordem,i,j,maux);
			}
			break; 
		case 1:
			j = num; 
			for(i=0; i<ordem; i++){
				MenorPrincipal (matriz,maux,i,j,ordem); 
				soma += matriz[i][j] * Cofator (ordem,i,j,maux);
			}
			break; 	
		}
	return soma; 
}
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
//----------------------------------------------------------
void SistemaTriangularInferior(int ordem, double matriz[][max], double b[], double x[]){
	double somatorio; 
	for(int i=0; i<ordem; i++){
		somatorio = 0; 
		for(int j=0; j<i; j++){
			somatorio += (x[j] * matriz[i][j]); 
		}
		x[i] = (b[i] - somatorio) / matriz[i][i];
	}
}
/*************************************************************
		     Cálculo dos Sistemas Genéricos
		     ------------------------------
					Decomposicao LU
*************************************************************/
void Uij(int linha, int coluna, double matriz[][max], double U[][max], double L[][max]){
	double somatorio=0;
	int k;
	for(k=0;k<=linha-1;k++){
		somatorio += (L[linha][k]*U[k][coluna]);
	}
	U[linha][coluna] = matriz[linha][coluna] - somatorio;
	return;
}
//-----------------------------------------------------------
void Lij(int linha, int coluna, double matriz[][max], double U[][max], double L[][max]){
	double somatorio=0;
	int k;
	for(k=0;k<=coluna-1;k++){
		somatorio += (L[linha][k]*U[k][coluna]);
	}
	L[linha][coluna]=(matriz[linha][coluna] - somatorio)/U[coluna][coluna];
	return;
}
//-----------------------------------------------------------
void ZeraMatriz(double matriz[max][max], int ordem){
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++){
			matriz[i][j]=0;
		}
	}
}
//-----------------------------------------------------------
void DecomposicaoLU(int ordem, double matriz[][max], double b[], double x[]){
	int i,j;
	double L[max][max], U[max][max],y[max]; //L = inferior ; U = superior
	ZeraMatriz(L, ordem);
	ZeraMatriz(U, ordem);
	for(i=0;i<ordem;i++){ //testar convergencia det(Ak)!= 0 ,  0 <= k < ordem
		if((Determinante(i+1,matriz)) == 0){
			printf("\nNao converge!");
			return;
		}
	}
	for(i=0;i<ordem;i++){ // variacao de linha
		for(j=0;j<ordem;j++){ //variacao de coluna
			if(i<=j){
				Uij(i,j,matriz,U,L);
				if(i==j){
					L[i][j]=1;
				}
			}
			else{
				Lij(i,j,matriz,U,L);
			}
		}
	}
	printf("\nMatriz U:");
	MostraMatriz(U,ordem);
	printf("\nMatriz L:");
	MostraMatriz(L,ordem);
	SistemaTriangularInferior(ordem,L,b,y);
	SistemaTriangularSuperior(ordem,U,y,x);
}
/*************************************************************
						Cholesky
*************************************************************/
int simetrica(double matriz[][max], int ordem){
	for(int i=0; i<ordem; i++){
		for(int j=i+1; j<ordem; j++){
			if(matriz[i][j]!=matriz[j][i]){
				return 0;
			}
		}
	}
	return 1;
}
//-----------------------------------------------------------
double L_Cholesky(double Aij, int i, int j, double L[][max]){
	double somatorio=0;
	
	if(i==j){ //elementos da diagonal
	 	for(int k=0; k<=i-1; k++){
			somatorio+=(pow(L[i][k], 2));
		}
		return sqrt(Aij-somatorio);
	}
	else{ //elementos fora da diagonal
		for(int k=0; k<=j-1; k++){
			somatorio+=( L[i][k] * L[j][k] );
		}
		return ((Aij-somatorio)/L[j][j]);
	}
}
//-----------------------------------------------------------
void Transpor_L(double L[][max], int ordem, double L_t[][max]){
	
	for(int i=0; i<ordem; i++){
		for(int j=i; j<ordem; j++){
			L_t[i][j]=L[j][i];
		}
	}
}
//-----------------------------------------------------------
void Cholesky(int ordem, double matriz[][max], double b[], double x[]){
	
	double L[max][max], L_t[max][max], y[max];
	ZeraMatriz(L, ordem);
	ZeraMatriz(L_t, ordem);
	//CONVERGENCIA:
	if(simetrica(matriz, ordem)==0){
		printf("O METODO NAO CONVERGE!");
	}
	for(int i=0; i<ordem; i++){
		if(Determinante(i+1, matriz)<=0){
			printf("\n O METODO NAO CONVERGE!");
			return;
		}
	}
	
	//METODO:
	
	for(int j=0; j<ordem; j++){
		for(int i=j; i<ordem; i++){
			L[i][j]=L_Cholesky(matriz[i][j], i, j, L);
			//printf("%lf \t", L[i][j]);
		}
		//printf("\n");
	}
	printf("\n Matriz L:\n");
	MostraMatriz(L, ordem);
	Transpor_L(L, ordem, L_t);
	printf("\n Matriz L_t:\n");
	MostraMatriz(L_t, ordem);
	system("pause");
	SistemaTriangularInferior(ordem, L, b, y);
	SistemaTriangularSuperior(ordem, L_t, y, x);
}
//-----------------------------------------------------------
/************************************************************
					Gauss Compacto
*************************************************************/
void TrocaLinhas(double matriz[max][max], int linha1, int linha2, int ordem){
	double aux;
	for(int j=0; j<ordem; j++){
		aux=matriz[linha1][j];
		matriz[linha1][j]=matriz[linha2][j];
		matriz[linha2][j]=aux;
	}
}
//-----------------------------------------------------------
void Uij_compacto(int linha, int coluna, double mExt[][max], double mLU[max][max]){
	double somatorio=0;
	int k;
	for(k=0;k<=linha-1;k++){
		somatorio += (mLU[linha][k]*mLU[k][coluna]);
	}
	mLU[linha][coluna] = mExt[linha][coluna] - somatorio;
}
//-----------------------------------------------------------
void Lij_compacto(int linha, int coluna, double mExt[][max], double mLU[max][max]){
	double somatorio=0;
	int k;
	for(k=0;k<=coluna-1;k++){
		somatorio += (mLU[linha][k]*mLU[k][coluna]);
	}
	mLU[linha][coluna]=(mExt[linha][coluna] - somatorio)/mLU[coluna][coluna];
}
//-----------------------------------------------------------
void GaussCompacto(int ordem, double matriz[][max], double b[], double x[]){
	double mExt[max][max], mLU[max][max];
	
	//CONVERGENCIA:
	for(int i=0; i<ordem; i++){
		if(Determinante(i+1, matriz)==0){
			printf("\n O METODO NAO CONVERGE!");
			return;
		}
	}
	//JUNTAR MATRIZES:
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem+1; j++){
			if(j<ordem){
				mExt[i][j]=matriz[i][j];
			}
			else{
				mExt[i][j]=b[i];
			}
		}
	}
	//METODO:
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem+1; j++){
			if(j>=i){
				Uij_compacto(i, j, mExt, mLU);
			}
			else{
				Lij_compacto(i, j, mExt, mLU);
			}
		}
	}
	for(int i=0; i<ordem; i++){
		b[i]=mLU[i][ordem];
	}
	SistemaTriangularSuperior(ordem, mLU, b, x);
}

/************************************************************
					    Gauss Jordan 
************************************************************/
void AmpliaMatriz(int ordem, double matriz[][max], double vetor[]){
	for(int i=0; i<ordem; i++)
		matriz[i][ordem] = vetor[i];
}
//-----------------------------------------------------------
void GaussJordan (int ordem, double matriz[][max], double b[], double x[]){
	
	int pivo;
	float m; 
	double maux[ordem][ordem+1]; 
	int coluna = 0; 
	
	if(Determinante(ordem,matriz) == 0){
		printf("\nMétodo não converge"); 
		system("pause"); 
		return; 
	}
	
	AmpliaMatriz(ordem,matriz,b); 
	printf("\n\n\tMatriz Ampliada\n"); 
	for(int i=0; i<ordem; i++){
		printf("\n\t"); 
		for(int j=0; j<=ordem; j++){
			printf(" %lf ",matriz[i][j]); 
		} 
	}
	
	for(coluna = 0; coluna < ordem; coluna++){
		
		//Escolha do pivo
		for(int i = coluna; i<ordem; i++){
			if(matriz[i][coluna] != 0){
				pivo = i; 
				break; 
			}
		}
		
		if(pivo != coluna){
			TrocaLinhas(matriz,pivo,coluna,ordem+1); 
		}
		
		for(int i=0; i<ordem; i++){
			if(i != coluna){				
				m = matriz[i][coluna] / matriz[coluna][coluna]; 	
				for(int j=0; j<ordem+1; j++){
					matriz[i][j] = matriz[i][j] - (m * matriz[coluna][j]); 
					if(j == ordem)
						b[i] = matriz[i][j]; //Atualiza o vetor b fora da matriz
				}
			}
		}	
	} 	
	
	SistemaTriangularSuperior(ordem,matriz,b,x);
	
	printf("\n\n\tMatriz Zerada\n"); 
	for(int i=0; i<ordem; i++){
		printf("\n\t"); 
		for(int j=0; j<=ordem; j++){
			printf(" %lf ",matriz[i][j]); 
		}
		
	} 
}

/*************************************************************
					     Matriz Inversa 
*************************************************************/
void GerarIdentidade (int ordem, double I[][max]){
	
	ZeraMatriz(I,ordem); 
	for(int i=0; i<ordem; i++)
		I[i][i] = 1; 	
}
//-----------------------------------------------------------
void CopiaMatriz(int ordem, double m1[][max], double m2[][max]){
	for(int i=0; i<ordem; i++)
		for(int j=0; j<ordem; j++)
			m2[i][j] = m1[i][j]; }
//-----------------------------------------------------------
void TransporMatriz(int ordem, double matriz[][max]){
	
	double maux[max][max]; 
	CopiaMatriz(ordem,matriz,maux); 
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++){
			matriz[i][j] = maux[j][i]; 
		}
	}
}
//-----------------------------------------------------------
void MatrizInversa(int ordem, double matriz[][max], double inversa[][max]){
	
	int op;
	double y[max][max]; 
	system("cls"); 
	printf("\n\tQual o método a ser utilizado?"); 
	printf("\n\t 1 - Decomposiçao LU"); 
	printf("\n\t 2 - Gauss Compacto"); 
	printf("\n\n\tOpção: "); 
	do{
		scanf("%d",&op); 
	}while(op != 1 && op != 2);
	
	GerarIdentidade(ordem,y); 	
	switch(op){
		case 1:				
			for(int i=0; i<ordem; i++)
				DecomposicaoLU(ordem,matriz,y[i],inversa[i]);					
			break; 
		case 2: 
			//for(int i=0; i<ordem; i++)
				//GaussCompacto(ordem,matriz,y[i],inversa[i]);
			break; 
	}	
	TransporMatriz(ordem,inversa); 
}
/*************************************************************
		    				Menu
**************************************************************/
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
//-----------------------------------------------------------
void ColetaDados(int *ordem, double matriz[][max], double b[max],  int opcao){
	
	// coleta dados para calculo do determinante
	if(opcao==1){ 
		//system("cls");
		printf("\n\tDigite a ordem da matriz:");
		scanf("%d", &(*ordem));
		printf("\nDigite a Matriz: \n"); 
		for(int i=0; i<(*ordem); i++){
			for(int j=0; j<(*ordem); j++){
				scanf("%lf",&matriz[i][j]);
			}
		}
	}	
	// coleta dados para Sistemas lineares
	if(opcao==2){
		system("cls");
		printf("\n\tDigite a ordem da matriz:");
		scanf("%d", &(*ordem));
		printf("\n\tDigite a Matriz: \n"); 
		for(int i=0; i<(*ordem); i++){
			for(int j=0; j<(*ordem); j++){
				scanf("%lf",&matriz[i][j]);
			}
		}
		printf("\n\tDigite o vetor de termos independentes:\n");
		for(int i=0; i<(*ordem); i++){
			scanf("%lf", &b[i]);
		}
	}
}
//-----------------------------------------------------------
void MostraMatriz(double matriz[max][max], int ordem){
	printf("\n\tMatriz:\n");
	for(int i=0; i<ordem; i++){
		for(int j=0; j<ordem; j++){
			printf("\t%.2lf", matriz[i][j]);
		}
		printf("\n");
	}
}
//-----------------------------------------------------------
void MostraSolucao(double x[max], int ordem){
	printf("\n\tSolucao:\n");
	for(int i=0; i<ordem; i++){
		printf("\t%lf\n", x[i]);
	}
}
//-----------------------------------------------------------
int main(){
	int opcao, ordem;
	double det, matriz[max][max], inv_matriz[max][max], b[max], solucao[max]; 
	
	/*freopen("matriz.txt","r",stdin); 	//Matriz para calcular gauss jordan
	//freopen("matriz2.txt","r",stdin); // Matriz superior
	//freopen("cholesky.txt", "r", stdin); // cholesky
	
	scanf("%d",&ordem); 
	for(int i=0; i<ordem; i++)
		for(int j=0; j<ordem; j++)
			scanf("%lf",&matriz[i][j]); 
	for(int i=0; i<ordem; i++)
		scanf("%lf",&b[i]); 
	
	MostraMatriz(matriz,ordem); 	
	GaussJordan(ordem,matriz,b,solucao); 
	//MostraMatriz()
	MostraSolucao(solucao,ordem); 
	system("pause"); 
	/**/	
	do{
		opcao = MenuMetodos();
		switch(opcao){
			case 1:		
				ColetaDados(&ordem, matriz, b, 1);
				MostraMatriz(matriz,ordem);
				det=Determinante(ordem, matriz); 
				printf("det(A) = %lf\n\t", det);
				system("pause"); 
				break;
			case 2: 
				ColetaDados(&ordem, matriz, b, 2);
				MostraMatriz(matriz, ordem);
				SistemaTriangularInferior(ordem, matriz, b, solucao); 
				MostraSolucao(solucao, ordem);
				system("pause");
				break; 
			case 3: 
				ColetaDados(&ordem, matriz, b, 2);
				MostraMatriz(matriz, ordem);
				SistemaTriangularSuperior(ordem, matriz, b, solucao); 
				MostraSolucao(solucao, ordem);
				system("pause");
				break; 
			case 4: 
				ColetaDados(&ordem, matriz, b, 2);
				MostraMatriz(matriz, ordem);
				DecomposicaoLU(ordem, matriz, b, solucao); 
				MostraSolucao(solucao, ordem);
				system("pause");
				break; 
			case 5: 
				ColetaDados(&ordem, matriz, b, 2);
				system("pause");
				MostraMatriz(matriz, ordem);
				Cholesky(ordem, matriz, b, solucao); 
				MostraSolucao(solucao, ordem);
				system("pause");
				break; 
			case 6:
				ColetaDados(&ordem, matriz, b, 2);
				MostraMatriz(matriz, ordem);
				GaussCompacto(ordem, matriz, b, solucao); 
				MostraSolucao(solucao, ordem);
				system("pause");
				break; 
			case 7: 
				ColetaDados(&ordem,matriz,b,2); 
				printf("\nDigite o vetor dos termos independentes: ");  
				GaussJordan(ordem,matriz,b,solucao); 
				MostraSolucao(solucao,ordem); 
				system("pause");  
				break; 
			case 8: 		//GPP sem troca de linhas 
				break; 
			case 9:			//Jacobi
				break; 
			case 10: 		//Gauss Seidel 
				break; 
			case 11: //Matriz inversa				
				ColetaDados(&ordem,matriz,b,1);
				MatrizInversa(ordem,matriz,inv_matriz);
				MostraMatriz(inv_matriz,ordem); 
				system("pause");  
		}
	}while(opcao != 12);
}
