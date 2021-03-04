#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    char url[101],carac;
    FILE *arq;
    int tam, max_it,opcao,i,j,k,t,n;//tam=tamanho,max_it=max iterações, n = tam-1, X[][]=matriz dos coeficientes, B[]=vetor dos termos independentes, X2[]= vetor dos resultados
    double erro_max, X[100][100],B[100],X2[100],X0[100],X1[100],m,erro,box,box2,soma1=0,soma2=0,aux,aux2,aux3;//tamanho 100, pois é o tamanho máximo possivel.

   do {
        printf("Escolha o metodo para o calculo da solucao\n");
        printf("1 - Sair\n2 - Calcular pelo Metodo de Gauss-Jordan com pivotacao parcial");
        printf("\n3 - Calcular pelo Metodo de Gauss-Jordan sem pivotacao\n4 - Metodo Iterativo de Jordan com Dominancia por linha\n5 - Metodo Iterativo de Seidel com Dominancia por linha\nOpcao:");
        scanf("%d", &opcao);
        if(opcao<1||opcao>5){
            printf("Opcao invalida.\nTente novamente!\n");
        }
    } while(opcao<1||opcao>5);
    //em caso de escolha 1 o programa fecha, em caso de qualquer outra opção ele continua rodando(o else abrange o programa inteiro).
    if(opcao == 1){
        printf("saindo do programa.");
        return 0;
    }
    else{
            fflush(stdin);
        printf("informe o nome do arquivo (com terminacao .txt e maximo de 100 caracteres):");
        scanf("%[^\n]s",url);
    arq = fopen(url,"r");//abrindo o arquivo com o nome dado pelo usuario.
    if(arq == NULL){
            printf("Erro, nao foi possivel abrir o arquivo\nReinicie o programa");
			return 0;
    }
    fscanf(arq,"%d %d %lf\n", &tam ,&max_it ,&erro_max);//escaneando a primeira linha para descobrir tamanho, max iterações e erro máximo.
     //função que escaneia o restante do arquivo e coloca seus valores em 2 matrizes (X é a matriz de coeficientes e B é a matriz de termos independentes).
     for(i=0 ;i<tam ;i++){
     for(j=0 ;j<tam ;j++){
    fscanf(arq,"%lf", &X[i][j]);
      }
    fscanf(arq,"%lf", &B[i]);
    }
fclose(arq);
     switch(opcao){
    case 2:
           n=tam-1;
           for(k=0;k<n;k++){
           for(i=k+1;i<=n;i++){
            if(fabs(X[i][k])>fabs(X[k][k])){
                for(t=0;t<=n+1;t++){
                    box=X[k][t];
                    X[k][t]=X[i][t];
                    X[i][t]=box;
                }
                    box2=B[k];
                    B[k]=B[i];
                    B[i]=box2;
            }
        }
        for(i=k+1;i<=n;i++){
            m=(X[i][k])/(X[k][k]);
            X[i][k]=0;
            for(j=k+1;j<=n;j++){
                X[i][j]=(X[i][j])-m*(X[k][j]);
            }
            B[i]=(B[i])-(m*(B[k]));
        }
    }

    for(k=n;k>0;k--){
        for(i=k-1;i>=0;i--){
            m=(X[i][k])/(X[k][k]);
            X[i][k]=0;
            for(j=k-1;j>=0;j--){
                X[i][j]=(X[i][j])-m*(X[k][j]);
            }
            B[i]=B[i]-(m*B[k]);
        }
    }



    for(i=0;i<=n;i++){
        X2[i]=(B[i])/(X[i][i]);
    }
        break;
    case 3:
        n=tam-1;
        for(k=0;k<n;k++){
        for(i=k+1;i<=n;i++){
            m=(X[i][k])/(X[k][k]);
            X[i][k]=0;
            for(j=k+1;j<=n;j++){
                X[i][j]=(X[i][j])-(m*(X[k][j]));
            }
            //b[i]=b[i]-m*b[k]
            B[i]=B[i]-(m*(B[k]));
            }
        }

        for(k=n;k>0;k--){
            for(i=k-1;i>=0;i--){
                m=(X[i][k])/(X[k][k]);
                X[i][k]=0;
                for(j=k-1;j>=0;j--){
                    X[i][j]=(X[i][j])-(m*(X[k][j]));
                }
                B[i]=(B[i])-(m*(B[k]));
            }
        }
        for(i=0;i<=n;i++){
            X2[i]=(B[i])/(X[i][i]);
        }
        break;
    case 4:
    for(k=0;k<tam;k++){//função para deixar a matriz dominante por linha
    j=-1;
    aux=X[k][k];
    for(i=0;i<tam;i++){
        if(fabs(aux)<fabs(X[k][i])){
            aux=X[k][i];
            j=i;
            }
        }
        if(fabs(aux)>fabs(X[k][k])){
                for(t=0;t<tam;t++){
                    aux2=X[k][t];
                    X[k][t]=X[j][t];
                    X[j][t]=aux2;
                }
            aux3=B[k];
            B[k]=B[j];
            B[j]=aux3;
            k--;
            }
     }
 for(i=0; i<tam; i++){
  X0[i]=1;
  }//colocando zero em todo o vetor, já que vamos considerar o vetor nulo como a aproximação inicial;
    printf("Considerando a aproximação inicial o vetor nulo.\n");
  k=0;
  do{
    for(i=0;i<tam;i++){
      soma1=0;
      for(j=0;j<tam;j++){
        if(j!=i){
        soma1=soma1+((X[i][j]*X0[j])/X[i][i]);
        }
      X1[i]=(B[i]/X[i][i])-soma1;
      }
    }

  //A partir daqui é a parte de erro máximo
    soma1=0;
    soma2=0;
  for(j=0;j<tam;j++){
    soma1=soma1+(((X1[j]-X0[j]))*(X1[j]-X0[j]));
    soma2=soma2+((X0[j])*(X0[j]));
  }
  soma1=sqrt(soma1);//calculado a norma dos vetores aqui
  soma2=sqrt(soma2);
  erro=soma1/soma2;
  //termina aqui
  for(i=0;i<tam;i++){
  X0[i]=X1[i];
  }
  k++;
  }while((k<max_it && k<=5000) && erro>=erro_max);
  for(i=0;i<tam;i++){
  X2[i]=X1[i];
  }
  break;
    case 5:

    for(k=0;k<tam;k++){//função para deixar a matriz dominante por linha.
    j=-1;
    aux=X[k][k];
    for(i=0;i<tam;i++){
        if(fabs(aux)<fabs(X[k][i])){
            aux=X[k][i];
            j=i;
            }
        }
        if(fabs(aux)>fabs(X[k][k])){
                for(t=0;t<tam;t++){
                    aux2=X[k][t];
                    X[k][t]=X[j][t];
                    X[j][t]=aux2;
                }
            aux3=B[k];
            B[k]=B[j];
            B[j]=aux3;
            k--;
            }
     }


for(i=0;i<tam;i++){//inicio do metodo
    X2[i]=0;
    X0[i]=0; //valor da aproximação inicial
 }
 k=0;
do{
 for(i=0;i<tam;i++){
   box=0;
   box2=0;
    for(j=0;j<i;j++){
     box=box+(X[i][j]*X2[j]);
   }
   for(j=(i+1);j<tam;j++){
    box2=box2+(X[i][j]*X0[j]);
   }
    X2[i]=(B[i]-box-box2)/X[i][i];
    //parte relacionada ao erro
    soma1=soma1+(((X2[i]-X0[i]))*(X2[i]-X0[i]));
    soma2=soma2+((X0[i])*(X0[i]));
  }
  for(i=0;i<tam;i++){
    X0[i]=X2[i];//atualizando o X0 para a proxima iteração
  }
  soma1=sqrt(soma1);//calculado a norma dos vetores aqui
  soma2=sqrt(soma2);
  erro=soma1/soma2;

 k++;
 }while((k<max_it && k<5000) && erro>erro_max);
        break;
        }
    }
    printf("\nO vetor solucao obtido e:\n");
    for(i=0;i<tam;i++){
    printf("%.4e\n", X2[i]);// %.4e é utilizado para imprimir o resultado em notação cientifica com 4 casas decimais
    }
    return 0;
}
//função para printar as matrizes para checar se o arquivo foi escaneado corretamente
 /* for(i=0;i<tam;i++){
    for(j=0;j<tam;j++){
        printf("%lf ",X[i][j]);
    }
    printf("%lf\n",B[i]);
}
    printf("\ntam:%d  max_it:%d  erro max:%lf\n",tam,max_it,erro_max);
*/
