//  Segundo Trabalho de MNC 
//  Sistemas Lineares
// Leonardo Oliveira

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <locale.h>


int menu(){
	
	int x; 
	
	printf("\nEscolha a opção: "); 
	printf("\n\n1 - blabla "); 
	printf("\n10 - Fechar programa"); 
	printf("\n\nOpção: "); 
	scanf("%d",&x); 
	
	return x; 
}

int main(){

	setlocale(LC_ALL,"Portuguese");
	int op; 
	
	do{
		op = menu(); 
		switch(op){
			case 1: 
					break; 
			case 2: 
					break; 
		}
	}while(op != 10); 
	
	printf("\nFim do programa"); 
}

