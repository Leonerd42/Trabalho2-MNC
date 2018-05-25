/***************************************************
		Trabalho de MNC - Sistemas Lineares
Lucas Henrique Russo do Nascimento, RA: 171025709
Bruna Lika Tamake, RA: 171024427
Leonardo Silva de Oliveira, RA: 171025903

*****************************************************/
#include<math.h>
#include<stdio.h>
#include<conio.h>
#include<ctype.h>
#include<stdlib.h>
#include <locale.h>

#define max 10

double Determinante (int ordem, double matriz[][max]);
void MostraMatriz(double matriz[max][max], int linha, int coluna);
void CopiaMatriz(int linha,int coluna, double m1[][max], double m2[][max]);

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
void SistemaTriangularSuperior(int ordem, double matriz[][max], double b[], double x[]){
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
int DecomposicaoLU(int ordem, double matriz[][max], double b[], double x[]){
	int i,j;
	double L[max][max], U[max][max],y[max]; //L = inferior ; U = superior
	ZeraMatriz(L, ordem);
	ZeraMatriz(U, ordem);
	for(i=0;i<ordem;i++){ //testar convergencia det(Ak)!= 0 ,  0 <= k < ordem
		if((Determinante(i+1,matriz)) == 0){
			printf("\nNao converge!");
			return 0;
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
	MostraMatriz(U,ordem,ordem);
	printf("\nMatriz L:");
	MostraMatriz(L,ordem,ordem);
	SistemaTriangularInferior(ordem,L,b,y);
	SistemaTriangularSuperior(ordem,U,y,x);
	return 1; 
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
int Cholesky(int ordem, double matriz[][max], double b[], double x[]){
	
	double L[max][max], L_t[max][max], y[max];
	ZeraMatriz(L, ordem);
	ZeraMatriz(L_t, ordem);
	//CONVERGENCIA:
	if(simetrica(matriz, ordem)==0){
		printf("\n\tO METODO NAO CONVERGE!");
		printf("\n\tMATRIZ NÃO É SIMETRICA"); 
		return 0; 
	}
	for(int i=0; i<ordem; i++){
		if(Determinante(i+1, matriz)<=0){
			printf("\n\tO METODO NAO CONVERGE!");
			printf("\n\tDETERMINANTE DE Ak = 0");
			return 0;
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
	MostraMatriz(L, ordem,ordem);
	Transpor_L(L, ordem, L_t);
	printf("\n Matriz L_t:\n");
	MostraMatriz(L_t, ordem,ordem);
	printf("\n\t"); 
	system("pause");
	SistemaTriangularInferior(ordem, L, b, y);
	SistemaTriangularSuperior(ordem, L_t, y, x);
	return 1; 
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
int GaussCompacto(int ordem, double matriz[][max], double b[], double x[]){
	double mExt[max][max], mLU[max][max];
	
	//CONVERGENCIA:
	for(int i=0; i<ordem; i++){
		if(Determinante(i+1, matriz)==0){
			printf("\n O METODO NAO CONVERGE!");
			return 0;
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
	return 1; 
}
/************************************************************
					    Gauss Jordan 
************************************************************/
void AmpliaMatriz(int ordem, double matriz[][max], double vetor[]){
	for(int i=0; i<ordem; i++)
		matriz[i][ordem] = vetor[i];
}
//-----------------------------------------------------------
int GaussJordan (int ordem, double matriz[][max], double b[], double x[]){
	
	int pivo;
	float m; 
	double maux[ordem][ordem+1]; 
	int coluna = 0; 
	
	if(Determinante(ordem,matriz) == 0){
		printf("\nMétodo não converge"); 
		system("pause"); 
		return 0; 
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
	
	printf("\n\n\tMatriz Escalonada\n"); 
	MostraMatriz(matriz,ordem,ordem+1); 
	return 1;
}
/***********************************************************
					GPP Sem troca de linhas
***********************************************************/
double Modulo (double numero){	
	return (numero >= 0) ? numero : numero*(-1); 	
}
//---------------------------------------------------------
int GaussPivoParcialSemTrocas(int ordem, double matriz[][max], double b[], double x[]){
	
	int coluna; 
	int posicao[ordem]; 
	int aux, aux2;  
	double m; 
	double maux[max][max]; 
	
	if(Determinante(ordem,matriz) == 0){
		printf("\nMétodo não converge\n"); 
		system("pause"); 
		return 0; 
	}
	
	AmpliaMatriz(ordem,matriz,b); 
	printf("\n\n\tMatriz Ampliada\n"); 
	for(int i=0; i<ordem; i++){
		printf("\n\t"); 
		for(int j=0; j<=ordem; j++)
			printf(" %.2lf ",matriz[i][j]); 
	}
	
	for(int i=0; i<ordem; i++)
		posicao[i] = i; 
	
	for(coluna=0; coluna<ordem; coluna++){
		
		//Escolha do pivo e troca das linhas no vetor posição 
		aux = coluna; 
		for(int i=coluna; i<ordem; i++)		
			if(Modulo(matriz[posicao[i]][coluna]) > Modulo(matriz[aux][coluna]))
				aux = i; 
		//Troca das posições do vetor posição 
		aux2 = posicao[coluna]; 
		posicao[coluna] = posicao[aux]; 
		posicao[aux] = aux2; 
		
		//Zerando a coluna atual 
		for(int i = coluna; i<ordem; i++){
	    	if(posicao[i] != posicao[coluna]){
		    	m = matriz[posicao[i]][coluna] / matriz[posicao[coluna]][coluna];
		    	for(int j = 0; j<ordem+1; j++){
		    		matriz[posicao[i]][j] -= m * matriz[posicao[coluna]][j];
				}
			}else continue;
		}
	}
	
	printf("\n\n\tMatriz escalonada\n");
	MostraMatriz(matriz,ordem,ordem+1); 
	
	//Copia matriz 
	CopiaMatriz(ordem,ordem+1,matriz,maux); 
	
	//Troca as linhas da matriz 
	for(int i = 0; i < ordem; i++){
		for(int j=0; j<ordem+1; j++){
			matriz[i][j] = maux[posicao[i]][j];
		}
	}	
	
	//Atualiza o vetor b na ordem correta
	for(int i=0; i<ordem; i++){
		b[i] = matriz[i][ordem]; 
	}
	SistemaTriangularSuperior(ordem,matriz,b,x);
	return 1;
}
/************************************************************
				  	   Métodos Iterativo 
				  	-----------------------
				  	       Jacobi
************************************************************/
int CriterioDasLinhas(int ordem, double matriz[max][max]){
	double maior=0, somatorio;
	for(int i=0; i<ordem; i++){
		somatorio=0;
		for(int j=0; j<ordem; j++){
			if(j==i){
				continue;
			}
			somatorio+=fabs(matriz[i][j]/matriz[i][i]);
		}
		if(i==0){
			maior=somatorio;
		}
		else{
			if(somatorio>maior){
				maior=somatorio;
			}
		}
	}
	if(maior>1){
		return 0;
	}
	return 1;
}
//----------------------------------------------------------
int CriterioDasColunas(int ordem, double matriz[max][max]){
	double maior=0, somatorio;
	for(int j=0; j<ordem; j++){
		somatorio=0;
		for(int i=0; i<ordem; i++){
			if(j==i){
				continue;
			}
			somatorio+=fabs(matriz[i][j]/matriz[j][j]);
		}
		if(j==0){
			maior=somatorio;
		}
		else{
			if(somatorio>maior){
				maior=somatorio;
			}
		}
	}
	if(maior>1){
		return 0;
	}
	return 1;
}
//----------------------------------------------------------
int Jacobi(int ordem, double matriz[max][max], double b[max], double x0[max], double E, int kmax, double x[max], int *k){
	double somatorio, maior_N, maior_D, norma;
	(*k) = 0; 
	
	//CONVERGENCIA:
	if(Determinante(ordem, matriz)==0){
		printf("\n\n\tMETODO NAO CONVERGE!\n");
		return 0;
	}
	for(int i=0; i<ordem; i++){
		if(matriz[i][i]==0){
			printf("\n\n\tMETODO NAO CONVERGE!\n");
			return 0;
		}
	}
	if(!CriterioDasLinhas(ordem, matriz)){
		if(!CriterioDasColunas(ordem, matriz)){
			printf("\n\n\tMETODO NAO CONVERGE!\n");
			return 0;
		}
	}
	
	//METODO:
	while((*k)<kmax){
		for(int i=0; i<ordem; i++){
			somatorio=0;
			for(int j=0; j<ordem; j++){
				if(i!=j){
					somatorio+=matriz[i][j]*x0[j];
				}
			}
			x[i]=(b[i]-somatorio)/matriz[i][i];
		}
		//CRITERIO DE PARADA:
		maior_N=fabs(x[0]-x0[0]);
		for(int i=1; i<ordem; i++){
			norma=fabs(x[i]-x0[i]);
			if(maior_N<norma){
				maior_N=norma;
			}
		}
		maior_D=fabs(x[0]);
		for(int i=1; i<ordem; i++){
			norma=fabs(x[i]);
			if(maior_D<norma){
				maior_D=norma;
			}
		}
		if((maior_N/maior_D)<E){
			return 1;
		}
		for(int i=0; i<ordem; i++){
			x0[i]=x[i];
		}
		(*k)++;
	}
	
	return 2; 
}
/************************************************************
					  	Gauss-Seidel
************************************************************/
int CriterioSassenfeld(int ordem, double matriz[][max]){
	int beta[max],i,j,k;
	double somatorio1, somatorio2, maior;
	for(i=0;i<ordem;i++){ //inicializando vetor
		beta[i]	=0;
	}
	for(i=0;i<ordem;i++){
		somatorio1=0;
		for(j=0;j<ordem;i++){
			somatorio1 += (fabs(matriz[i][j]/matriz[i][i])*beta[j]);
		}
		for(k=i+1;k<ordem;i++){
			somatorio2 += matriz[i][k]/matriz[i][i];
		}
		beta[i] = somatorio1 + somatorio2;
	}
	maior = beta[0];
	for(i=1; i<ordem; i++){
		if(beta[i]> maior)
			maior = beta[i];
	}
	if(maior>1){ //nao converge
		return 0;
	}
	else //converge
		return 1;
}
//-----------------------------------------------------------
int GaussSeidel(int ordem, double matriz[][max], double b[max], double x0[],float E, int kmax, double x[max], int *k){
	double somatorio, maior_N, maior_D, norma;
	
	(*k) = 0; 
	
	if(Determinante(ordem, matriz)==0){ //det(A) !=0
		printf("\n\n\tMETODO NAO CONVERGE!\n");
		return 0;
	}
	for(int i=0; i<ordem; i++){
		if(matriz[i][i]==0){ // matriz[i][i] !=0
			printf("\n\n\tMETODO NAO CONVERGE!\n");
			return 0;
		}
	}
	if(!CriterioDasLinhas(ordem, matriz)){
		if(!CriterioSassenfeld(ordem, matriz)){
			printf("\n\n\tMETODO NAO CONVERGE!\n");
			return 0;
		}
	}
	while((*k)<kmax){
		for(int i=0; i<ordem; i++){
			somatorio=0;
			for(int j=0; j<ordem; j++){
				if(i<j){
					somatorio+=matriz[i][j]*x0[j];
				}
				else{
					if(i>j){
						somatorio+=matriz[i][j]*x[j];
					}
				}
			}
			x[i]=(b[i]-somatorio)/matriz[i][i];	
		}
		maior_N=fabs(x[0]-x0[0]);
		for(int i=1; i<ordem; i++){
			norma=fabs(x[i]-x0[i]); //norma de x_atual - x_antigo
			if(maior_N<norma){
				maior_N=norma;
			}
		}
		maior_D=fabs(x[0]);
		for(int i=1; i<ordem; i++){
			norma=fabs(x[i]); //norma de x_atual
			if(maior_D<norma){
				maior_D=norma;
			}
		}
		if((maior_N/maior_D)<E){ //criterio de parada
			return 1;
		}
		for(int i=0; i<ordem; i++){
			x0[i]=x[i];
		}
		(*k)++;
	}
	return 2; 
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
void CopiaMatriz(int linha,int coluna, double m1[][max], double m2[][max]){
	for(int i=0; i<linha; i++)
		for(int j=0; j<coluna; j++)
			m2[i][j] = m1[i][j]; }
//-----------------------------------------------------------
void TransporMatriz(int ordem, double matriz[][max]){
	
	double maux[max][max]; 
	CopiaMatriz(ordem,ordem,matriz,maux); 
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
			for(int i=0; i<ordem; i++)
				GaussCompacto(ordem,matriz,y[i],inversa[i]);
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
		printf("\n\tMENU - Sempre que for digitar um numero racional, use a virgula\n");
		//printf("\n\t1-Determinante \n\t2-Sistema Triangular Inferior \n\t3-Sistema Triangular Superior \n\t4-Decomposicao LU \n\t5-Cholesky \n\t6-Gauss Compacto \n\t7-Gauss Jordan \n\t8-GPP SEM troca de linhas \n\t9-Jacobi \n\t10-Gauss Seidel \n\t11-Matriz Inversa\n\t12-Fechar programa\n");
		printf("\n\t 1 - Decomposição LU"); 
		printf("\n\t 2 - Cholesky"); 
		printf("\n\t 3 - Gauss Compacto "); 
		printf("\n\t 4 - Gauss Jordan"); 
		printf("\n\t 5 - GPP Sem troca de Linhas"); 
		printf("\n\t 6 - Jacobi"); 
		printf("\n\t 7 - Gauss Seidel "); 
		printf("\n\t 8 - Matriz Inversa"); 
		printf("\n\t 9 - Fechar Programa"); 
		printf("\n\n\tDigite a opcao:");
		scanf("%d",&op);
		if(op<1 || op>9){
			printf("\n\nOpcao Invalida!\n");
		}
	}while(op<1 || op>9);
	return op;
}
//-----------------------------------------------------------
void ColetaDados(int *ordem, double matriz[][max], double b[max],  int opcao){
	
	// coleta dados para calculo do determinante
	switch(opcao){
		case 1:
			//system("cls");
			printf("\n\tDigite a ordem da matriz: ");
			scanf("%d", &(*ordem));
			printf("\nDigite a Matriz: \n"); 
			for(int i=0; i<(*ordem); i++)
				for(int j=0; j<(*ordem); j++)
					scanf("%lf",&matriz[i][j]);				
			break; 
		case 2: 
			//system("cls");
			printf("\n\tDigite a ordem da matriz: ");
			scanf("%d", &(*ordem));
			printf("\n\tDigite a Matriz: \n"); 
			for(int i=0; i<(*ordem); i++)
				for(int j=0; j<(*ordem); j++)
					scanf("%lf",&matriz[i][j]);				
			printf("\n\tDigite o vetor de termos independentes:\n");
			for(int i=0; i<(*ordem); i++)
				scanf("%lf", &b[i]);
			break; 
	}
}
//-----------------------------------------------------------
void ColetaDadosIterativos(int *ordem, double matriz[][max], double b[], double x0[], float *e, int *max_ite){
	
	//system("cls"); 
	printf("\n\tDigite a ordem da matriz: "); 
	scanf("%d",&(*ordem)); 
	printf("\n\tDigite a matriz: \n"); 
	for(int i=0; i<(*ordem); i++)
		for(int j=0; j<(*ordem); j++)
			scanf("%lf",&matriz[i][j]); 
	printf("\n\tDigite o vetor dos termos independentes:\n"); 
	for(int i=0; i<(*ordem); i++)
		scanf("%lf",&b[i]); 
	printf("\n\tDigite o vetor da aproximação inicial:\n"); 
	for(int i=0; i<(*ordem); i++)
		scanf("%lf",&x0[i]); 
	printf("\n\tDigite a precisão desejada: "); 
	scanf("%f",&(*e)); 
	printf("\n\tDigite o número maximo de iterações: ");
	scanf("%d",&(*max_ite)); 
	
}
//-----------------------------------------------------------
void MostraMatriz(double matriz[max][max], int linha, int coluna){
	printf("\n");
	for(int i=0; i<linha; i++){
		for(int j=0; j<coluna; j++){
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
	setlocale(LC_ALL,"Portuguese");
	int opcao, ordem, max_ite, iteracoes; 
	int convergencia; 
	float precisao; 
	double det, matriz[max][max], inv_matriz[max][max];	//Matrizes
	double x0[max], b[max], solucao[max]; //Vetores
		
	do{
		opcao = MenuMetodos();
		switch(opcao){
			case 1: 			//Decomposição LU
				system("cls"); 
				printf("\n\n\tMétodo da Decomposição LU\n\n"); 
				ColetaDados(&ordem, matriz, b, 2);
				MostraMatriz(matriz, ordem,ordem);
				convergencia = DecomposicaoLU(ordem, matriz, b, solucao); 
				if(convergencia)
					MostraSolucao(solucao, ordem);
				printf("\n\t");
				system("pause");
				break; 
			case 2: 			//Cholesky
				system("cls"); 
				printf("\n\n\tMétodo de Cholesky\n\n"); 
				ColetaDados(&ordem, matriz, b, 2);
				MostraMatriz(matriz, ordem,ordem);
				convergencia = Cholesky(ordem, matriz, b, solucao); 
				if(convergencia)
					MostraSolucao(solucao, ordem);
				printf("\n\t");
				system("pause");
				break; 
			case 3: 			//Gauss compacto
				system("cls"); 
				printf("\n\n\tMétodo Gauss Compacto\n\n");
				ColetaDados(&ordem, matriz, b, 2);
				MostraMatriz(matriz, ordem,ordem);
				convergencia = GaussCompacto(ordem, matriz, b, solucao); 
				if(convergencia)
					MostraSolucao(solucao, ordem);
				printf("\n\t");
				system("pause");
				break; 
			case 4: 			//Gauss Jordan
				system("cls"); 
				printf("\n\n\tMétodo Gauss Jordan \n\n");
				ColetaDados(&ordem,matriz,b,2);  
				convergencia = GaussJordan(ordem,matriz,b,solucao); 
				if(convergencia)
					MostraSolucao(solucao,ordem); 
				printf("\n\t");
				system("pause");  
				break; 
			case 5: 			//GPP Sem troca de linhas 
				system("cls"); 
				printf("\n\n\tMétodo de Gauss Pivo Parcial Sem Troca de Linhas \n\n");
				ColetaDados(&ordem,matriz,b,2); 
				convergencia = GaussPivoParcialSemTrocas(ordem,matriz,b,solucao); 
				if(convergencia)
					MostraSolucao(solucao,ordem); 
				printf("\n\t");
				system("pause"); 
				break; 
			case 6:				//Jacobi
				system("cls"); 
				printf("\n\n\tMétodo Iterativo - Jacobi\n\n");
				ColetaDadosIterativos(&ordem,matriz,b,x0,&precisao,&max_ite); 
				convergencia = Jacobi(ordem,matriz,b,x0,precisao,max_ite,solucao,&iteracoes); 
				if(convergencia == 1){
					MostraSolucao(solucao,ordem); 
					printf("\n\tNúmero de iterações: %d\n",iteracoes); 
				}else if(convergencia == 2){
					printf("\n\tNúmero máximo de iterações atingidas!"); 
					printf("\n\tResultado obtido dentro das iterações:\n"); 
					MostraSolucao(solucao,ordem); 
				}
				printf("\n\t");
				system("pause"); 
				break; 
			case 7: 			//Gauss Seidel 
				system("cls"); 
				printf("\n\n\tMétodo Iterativo - Gauss Seidel \n\n");
				ColetaDadosIterativos(&ordem,matriz,b,x0,&precisao,&max_ite); 
				convergencia = GaussSeidel(ordem,matriz,b,x0,precisao,max_ite,solucao,&iteracoes); 
				if(convergencia == 1){
					MostraSolucao(solucao,ordem); 
					printf("\nNúmero de iterações: %d\n",iteracoes); 
				}else if(convergencia == 2){
					printf("\n\tNúmero máximo de iterações atingidas!"); 
					printf("\n\tResultado obtido dentro das iterações:\n"); 
					MostraSolucao(solucao,ordem); 
				}
				printf("\n\t");
				system("pause"); 
				break; 
			case 8: 			//Matriz inversa	
				system("cls"); 
				printf("\n\n\tCálculo da matriz inversa\n\n");			
				ColetaDados(&ordem,matriz,b,1);
				MatrizInversa(ordem,matriz,inv_matriz);
				MostraMatriz(inv_matriz,ordem,ordem); 
				printf("\n\t");
				system("pause");  
		}
	}while(opcao != 9);
	
	printf("\n\tFim do programa!\n\n"); 
	printf("\n\tPrograma Desenvolvido por:"); 
	printf("\n\t Bruna Lika Tamake");
	printf("\n\t Leonardo Silva de Oliveira");
	printf("\n\t Lucas Henrique Russo do Nascimento");
	printf("\n\n\t3° Termo - BCC - UNESP BAURU\n\n");
}
