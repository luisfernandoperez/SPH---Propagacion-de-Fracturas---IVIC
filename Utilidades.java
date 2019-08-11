import java.io.*;
import java.util.*;
import java.lang.*;
import java.util.stream.*;
import java.util.Map.Entry;
public class Utilidades
{
	Parametros Param 	= new Parametros();
  
    public double Producto_Punto(double[] vectorA, double[] vectorB)
    {
		int longitud_a = vectorA.length;
		int longitud_b = vectorB.length;
		if(longitud_a != longitud_b  || longitud_a < 2 || longitud_b < 2)
		{
			System.out.println("La dimensión de los vectores son diferentes dentre si y no se puede computar el producto punto o las dimensiones de los mismos son menores a 2");
			Limpiar_Archivos(new File(Param.directorio_base), "class");
			System.exit(1);
			return 1.0;
		}
		else
		{
			if(longitud_a == 2)
				return vectorA[0]*vectorB[0] + vectorA[1]*vectorB[1];
			else
				return vectorA[0]*vectorB[0] + vectorA[1]*vectorB[1] + vectorA[2]*vectorB[2];
		}
    }
    
    public boolean Cambio(boolean a, boolean b)
    {
		return (a || b);
	}
    
    
    public double[] Producto_Cruz(double[] vectorA, double[] vectorB)
    {
		int longitud_a = vectorA.length;
		int longitud_b = vectorB.length;
		double[] resultado = new double[longitud_a];
		if(longitud_a != longitud_b  || longitud_a < 2 || longitud_b < 2)
		{
			System.out.println("La dimensión de los vectores son diferentes dentre si y no se puede computar el producto cruz o las dimensiones de los mismos son menores a 2");
			Limpiar_Archivos(new File(Param.directorio_base), "class");
			System.exit(1);
			return resultado;
		}
		else
		{
			if(longitud_a == 2)
			{
				resultado = new double[longitud_a - 1];
				resultado[0] = vectorA[0]*vectorB[1] - vectorA[1]*vectorB[0];
			}
			if(longitud_a == 3)
			{	
				resultado[0] = vectorA[1]*vectorB[2] - vectorA[2]*vectorB[1];
				resultado[1] = vectorA[2]*vectorB[0] - vectorA[0]*vectorB[2];
				resultado[2] = vectorA[0]*vectorB[1] - vectorA[1]*vectorB[0];
			}
			
		}
		return resultado;
    }
    
    public double[] Producto_Matriz_Vector(double[][] matriz, double[] vector)
    {
        double[] Resultado = new double[matriz.length];
        if(matriz[0].length == vector.length)
        {
            for(int i = 0; i< matriz.length; i++)
            {
                for(int j = 0; j< vector.length; j++)
                {
                	Resultado[i] += matriz[i][j]*vector[j];
                }
            }
        }
        else
        {
           System.out.println("La dimensión del vector es diferente del número de columnas de la matriz");
        }
        
        return Resultado;
    }

	public double[][] Producto_Vectores(double[] vectorA, double[] vectorB)
	{
		double[][] Resultado = new double[vectorA.length][vectorB.length];
		for(int i =0; i< vectorA.length; i++)
		{
			for(int j =0; j< vectorB.length; j++)
			{
				Resultado[i][j] = vectorA[i]*vectorB[j];
			}
		}
		return Resultado;
	}
    
	public double[][] Producto_Matrices(double[][] matrizA, double[][] matrizB)
	{
		double[][] Resultado = new double[matrizA.length][matrizB[0].length];
		if(matrizA[0].length == matrizB.length)
		{
			for(int i =0; i< matrizA.length; i++)
			{
				for(int j =0; j< matrizB[0].length; j++)
				{
					for(int k =0; k< matrizB.length; k++)
					{
						Resultado[i][j] += matrizA[i][k]*matrizB[k][j];
					}
				}
			}
		}
		else
		{
			System.out.println("La columna de la primera matriz es diferente a la fila de la segunda");
		}
		return Resultado;
	}

	public double Doble_Producto_Punto_Matrices(double[][] matrizA, double[][] matrizB)
	{
		double Resultado = 0;
		if(matrizA.length == matrizB.length && matrizA[0].length == matrizB[0].length)
		{
			for(int i =0; i< matrizA.length; i++)
			{
				for(int j =0; j< matrizA[0].length; j++)
				{
					Resultado += matrizA[i][j]*matrizB[i][j];
				}
			}
		}
		else
		{
			System.out.println("Las matrices no tienen la misma dimension, metodo Doble_Producto_Punto_Matrices");
		}
		return Resultado;
	}

	public double[] Producto_Vector_Escalar(double escalar, double[] vector)
    {
        double[] Resultado = new double[vector.length];
        for(int i =0; i< vector.length; i++)
		{
			Resultado[i] = escalar*vector[i];
		}
        return Resultado;
    }

	public double[][] Producto_Matriz_Escalar(double escalar, double[][] Matriz)
    {
        double[][] Resultado = new double[Matriz.length][Matriz[0].length];
        for(int i =0; i< Matriz.length; i++)
			{
				for(int j =0; j< Matriz[0].length; j++)
				{
					Resultado[i][j] = escalar*Matriz[i][j];
				}
			}
        return Resultado;
    }

	public double[][] Suma_Matrices(double[][] matrizA, double[][] matrizB)
	{
		double[][] Resultado = new double[matrizA.length][matrizA[0].length];
		if(matrizA.length == matrizB.length && matrizA[0].length == matrizB[0].length)
		{
			for(int i = 0; i < matrizA.length; i++)
			{
				for(int j = 0; j < matrizA[0].length; j++)
				{
					Resultado[i][j] = matrizA[i][j] + matrizB[i][j];
				}
			}
		}
		else
		{
			System.out.println("Las matrices no tienen la misma dimension, metodo Suma_Matrices");
		}
		return Resultado;
	}


	public double[][] Resta_Matrices(double[][] matrizA, double[][] matrizB)
	{
		double[][] Resultado = new double[matrizA.length][matrizA[0].length];
		if(matrizA.length == matrizB.length && matrizA[0].length == matrizB[0].length)
		{
			for(int i = 0; i < matrizA.length; i++)
			{
				for(int j = 0; j < matrizA[0].length; j++)
				{
					Resultado[i][j] = matrizA[i][j] - matrizB[i][j];
				}
			}
		}
		else
		{
			System.out.println("Las matrices no tienen la misma dimension, metodo Resta_Matrices");
		}
		return Resultado;
	}
	
	public int[][] Resta_Matrices(int[][] matrizA, int[][] matrizB)
	{
		int[][] Resultado = new int[matrizA.length][matrizA[0].length];
		if(matrizA.length == matrizB.length && matrizA[0].length == matrizB[0].length)
		{
			for(int i = 0; i < matrizA.length; i++)
			{
				for(int j = 0; j < matrizA[0].length; j++)
				{
					Resultado[i][j] = matrizA[i][j] - matrizB[i][j];
				}
			}
		}
		else
		{
			System.out.println("Las matrices no tienen la misma dimension, metodo Resta_Matrices");
		}
		return Resultado;
	}

	public double[] Suma_Vectores(double[] vectorA, double[] vectorB)
	{
		double[] Resultado = new double[vectorA.length];
		if(vectorA.length == vectorB.length)
		{
			for(int i =0; i< vectorA.length; i++)
			{
				Resultado[i] = vectorA[i] + vectorB[i];
			}
		}
		else
		{
			System.out.println("Los vectores no tienen la misma dimension");
		}
		return Resultado;
	}

    public int Suma_Componente_Vector(int[] vectorA)
    {
        int Resultado = 0;
        for(int i =0; i< vectorA.length; i++)
        {
            Resultado += vectorA[i];
        }
        return Resultado;
    }

    public double Suma_Componente_Vector(double[] vectorA)
    {
        double Resultado = 0.0;
        for(int i =0; i< vectorA.length; i++)
        {
            Resultado += vectorA[i];
        }
        return Resultado;
    }

    public int[] Suma_Segunda_Componente_Matriz(int[][] vectorA)
    {
        int[] Resultado = new int[vectorA.length];
        for(int i =0; i< vectorA.length; i++)
        {
			for(int j =0; j< vectorA[i].length; j++)
			{
				Resultado[i] += vectorA[i][j];
			}
        }
        return Resultado;
    }

    public double[] Suma_Segunda_Componente_Matriz(double[][] vectorA)
    {
        double[] Resultado = new double[vectorA.length];
        for(int i =0; i< vectorA.length; i++)
        {
			for(int j =0; j< vectorA[i].length; j++)
			{
				Resultado[i] += vectorA[i][j];
			}
        }
        return Resultado;
    }
    
    public double[] Suma_Primera_Componente_Matriz(double[][] vectorA)
    {
        double[] Resultado = new double[vectorA[0].length];
        for(int i =0; i< vectorA[0].length; i++)
        {
			for(int j =0; j< vectorA.length; j++)
			{
				Resultado[i] += vectorA[j][i];
			}
        }
        return Resultado;
    }
       
    public double[] Union_Tensor_rango_1(double[] vectorA, double[] vectorB)
    {
        double[] union_vectores;

        int k =0, q=0;
        int largo_efectivo = vectorA.length + vectorB.length;
        union_vectores = new double[largo_efectivo];

        if(vectorA.length != 0)
        {
            for(int i =0; i< vectorA.length; i++)
            {
                union_vectores[i] = vectorA[i];
            }
        }

        for(int j = vectorA.length; j< largo_efectivo; j++)
        {
            union_vectores[j] = vectorB[j - vectorA.length];
        }

        return union_vectores ;
    }

    public int[] Union_Array(int[] vectorA, int[] vectorB)
    {
        int[] union_vectores;

        int k =0, q=0;
        int largo_efectivo = vectorA.length + vectorB.length;
        union_vectores = new int[largo_efectivo];

        if(vectorA.length != 0)
        {
            for(int i =0; i< vectorA.length; i++)
            {
                union_vectores[i] = vectorA[i];
            }
        }

        for(int j = vectorA.length; j< largo_efectivo; j++)
        {
            union_vectores[j] = vectorB[j - vectorA.length];
        }

        return union_vectores ;
    }
    
    public double[][] Union_Tensor_rango_2(double[][] MatrizA, double[][] MatrizB)
    {
        int ancho = 1;
        if(MatrizA.length == 0)
        {
            ancho = Parametros.dimension;
        }
        else
        {
            ancho = MatrizA[0].length;
        }
        double[][] union_matrices;
        if(ancho == MatrizB[0].length)
        {
            int largo_efectivo = MatrizA.length + MatrizB.length;
            union_matrices = new double[largo_efectivo][ancho];

            if(MatrizA.length != 0)
            {
                for(int i =0; i< MatrizA.length; i++)
                {
                    for(int k =0; k< ancho; k++)
                    {
                        union_matrices[i][k] = MatrizA[i][k];
                    }
                }
            }

            for(int j = MatrizA.length; j< largo_efectivo; j++)
            {
                for(int k =0; k< ancho; k++)
                {
                    union_matrices[j][k] = MatrizB[j - MatrizA.length][k];
                }
            }
        }
        else
        {
            union_matrices = new double[0][0];
            System.out.println("Los dos objetos no tienen la misma dimension");
        }
        return union_matrices;
    }

    public double[][][] Union_Tensor_rango_3(double[][][] MatrizA, double[][][] MatrizB)
    {
        int ancho =  1;
        int profundidad = 1;
        if(MatrizA.length == 0)
        {
            ancho = Parametros.dimension;
            profundidad = ancho;
        }
        else
        {
            ancho = MatrizA[0].length;
            profundidad = MatrizA[0][0].length;
        }
        double[][][] union_matrices;
        if(ancho == MatrizB[0].length && profundidad ==  MatrizB[0][0].length)
        {
            int largo_efectivo = MatrizA.length + MatrizB.length;
            union_matrices = new double[largo_efectivo][ancho][profundidad];

            if(MatrizA.length != 0)
            {
                for(int i =0; i< MatrizA.length; i++)
                {
                    for(int k =0; k< ancho; k++)
                    {
                        for(int u =0; u< profundidad; u++)
                        {
                            union_matrices[i][k][u] = MatrizA[i][k][u];
                        }
                    }
                }
            }

            for(int j = MatrizA.length; j< largo_efectivo; j++)
            {
                for(int k =0; k< ancho; k++)
                {
                    for(int u =0; u< profundidad; u++)
                    {
                        union_matrices[j][k][u] = MatrizB[j - MatrizA.length][k][u];
                    }
                }
            }
        }
        else
        {
            union_matrices = new double[0][0][0];
            System.out.println("Los dos objetos no tienen la misma dimension");
        }
        return union_matrices;
    }


	public double[] Resta_Vectores(double[] vectorA, double[] vectorB)
	{
		double[] Resultado = new double[vectorA.length];
		if(vectorA.length == vectorB.length)
		{
			for(int i =0; i< vectorA.length; i++)
			{
				Resultado[i] = vectorA[i] - vectorB[i];
			}
		}
		else
		{
			System.out.println("Los vectores no tienen la misma dimension");
		}
		return Resultado;
	}

    public double[][] Inversa_Matriz_2x2_3x3(int k, double[][] matriz)
    {
        int dimension 		= matriz.length;
        double[][] inversa 	= new double[dimension][matriz[0].length];
        double determinante = Determinante_R2_2x2_3x3(matriz);
        
        if(dimension != matriz[0].length)
        {
            System.out.println("La matriz no es cuadrada y no se puede invertir");
            System.exit(1);
        }
        
        else
        {
            if(dimension == 2)
            {            
                if(determinante != 0.)
                {
                    inversa[0][0] =   matriz[1][1]/determinante;
                    inversa[0][1] = - matriz[0][1]/determinante;
                    inversa[1][0] = - matriz[1][0]/determinante;
                    inversa[1][1] =   matriz[0][0]/determinante;
                }
                
                else
                {
					inversa[0][0] = 1;
                    inversa[0][1] = 0;
                    inversa[1][0] = 0;
                    inversa[1][1] = 1;
                    System.out.println("Determinante nulo, no se puede invertir la matriz 2x2 de la particula " + k);
                     //System.exit(1);
                }
                
            }
            
            if(dimension == 3)
            {
                
                if(determinante != 0)
                {
                    inversa[0][0] = (matriz[1][1]*matriz[2][2] - matriz[1][2]*matriz[2][1])/determinante;
                    inversa[0][1] = (matriz[0][2]*matriz[2][1] - matriz[0][1]*matriz[2][2])/determinante;
                    inversa[0][2] = (matriz[0][1]*matriz[1][2] - matriz[0][2]*matriz[1][1])/determinante;
                    inversa[1][0] = (matriz[1][2]*matriz[2][0] - matriz[1][0]*matriz[2][2])/determinante;
                    inversa[1][1] = (matriz[0][0]*matriz[2][2] - matriz[0][2]*matriz[2][0])/determinante;
                    inversa[1][2] = (matriz[0][2]*matriz[1][0] - matriz[0][0]*matriz[1][2])/determinante;
                    inversa[2][0] = (matriz[1][0]*matriz[2][1] - matriz[1][1]*matriz[2][0])/determinante;
                    inversa[2][1] = (matriz[0][1]*matriz[2][0] - matriz[0][0]*matriz[2][1])/determinante;
                    inversa[2][2] = (matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0])/determinante;
                }
                
                else
                {
				for(int i =0;  i<dimension; i++)
					for(int j =0;  j<dimension; j++)
						inversa[i][j] = Delta_Kronecker(i,j);
                    //System.out.println("Determinante nulo, no se puede invertir la matriz 3x3 de la particula " + k);
                    //Limpiar_Archivos(new File(Param.directorio_base), "class");
                    //System.exit(1);
                }
            }
            
            if(dimension != 2 && dimension != 3)
            {
                System.out.println("Este metodo solo calcula la inversa de matrices 2x2 y 3x3");
            }
        }
        return inversa;
    }
    
    public double Traza_R2(double[][] matriz)
    {
		double traza = 0;
		if(matriz.length == matriz[0].length)
		{
			for(int i=0; i < matriz.length; i++) traza += matriz[i][i];
		}
        else System.out.println("No se puede sacar traza a una matriz que no sea cuadrada");
        return traza;
    }
    
    public double Determinante_R2_2x2_3x3(double[][] matriz)
    {
		int dimension = matriz.length;
        double determinante = 0;
        
        if(dimension != matriz[0].length)
        {
            System.out.println("La matriz no es cuadrada y no se puede obtener determinante");
        }
        if(dimension == 2)
        {
			determinante =  matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0];
		}
		if(dimension == 3)
        {
			determinante = matriz[0][0]*matriz[1][1]*matriz[2][2] + matriz[0][1]*matriz[1][2]*matriz[2][0] + matriz[0][2]*matriz[1][0]*matriz[2][1] - (matriz[2][0]*matriz[1][1]*matriz[0][2] + matriz[2][1]*matriz[1][2]*matriz[0][0] + matriz[2][2]*matriz[1][0]*matriz[0][1]);
		}
		if(dimension != 2 && dimension != 3)
        {
			 System.out.println("Este metodo solo calcula la inversa de matrices 2x2 y 3x3");
		}
		return determinante;
	}

	public double[] Vector_de_Traza_R2(double[][][] matriz)
    {
		int dimension = matriz[0].length;
		double[] traza = new double[matriz.length];
		if(dimension == matriz[0][0].length)
		{
			for(int j=0; j < matriz.length; j++) 
			{
				for(int i=0; i < dimension; i++) traza[j] += matriz[j][i][i];
			}
		}
        else System.out.println("No se puede sacar traza a una matriz que no sea cuadrada");
        return traza;
    }
    
    public double[][] Transpuesta(double[][] matriz)
    {
        double[][] transpuesta = new double[matriz[0].length][matriz.length];
        
        for(int i= 0; i<matriz[0].length; i++)
        {
            for(int j= 0; j<matriz.length; j++)
            {
                transpuesta[i][j] = matriz[j][i];
            }
        }
        return transpuesta;
    }
    
    public int[][] Transpuesta(int[][] matriz)
    {
        int[][] transpuesta = new int[matriz[0].length][matriz.length];
        
        for(int i= 0; i<matriz[0].length; i++)
        {
            for(int j= 0; j<matriz.length; j++)
            {
                transpuesta[i][j] = matriz[j][i];
            }
        }
        return transpuesta;
    }
    
    
       public double[][][] Transpuesta_Primer_Indice_Al_Final(double[][][] tensor)
    {
        double[][][] transpuesta = new double[tensor[0].length][tensor[0][0].length][tensor.length];
        
        for(int i= 0; i<tensor[0].length; i++)
        {
			 for(int j= 0; j<tensor[0][0].length; j++)
			{
				for(int k= 0; k<tensor.length; k++)
				{
					transpuesta[i][j][k] = tensor[k][i][j];
				}
			}
        }
        return transpuesta;
    }

	public double[][] Matriz_Identidad_Rango_2(int dimension)
    {
        double[][] Resultado = new double[dimension][dimension];
        for(int i =0; i< dimension; i++)
		{
			for(int j =0; j< dimension; j++)
			{
				Resultado[i][j] = Delta_Kronecker(i, j);
			}
		}
        return Resultado;
    }

	public double[][] Tensor_De_Unos_Rango_2(int dimension)
    {
        double[][] Resultado = new double[dimension][dimension];
        for(int i =0; i< dimension; i++)
		{
			for(int j =0; j< dimension; j++)
			{
				Resultado[i][j] = 1.0;
			}
		}
        return Resultado;
    }
    
    public double[] Tensor_Cero_Rango_1(int a)
    {
        double[] Resultado = new double[a];
        for(int i =0; i< a; i++)
		{
			Resultado[i] = 0;
		}
        return Resultado;
    }

	public double[][] Tensor_Cero_Rango_2(int a, int b)
    {
        double[][] Resultado = new double[a][b];
        for(int i =0; i< a; i++)
		{
			for(int j =0; j< b; j++)
			{
				Resultado[i][j] = 0;
			}
		}
        return Resultado;
    }
    
    	public double[][][] Tensor_Cero_Rango_3(int a, int b, int c)
    {
        double[][][] Resultado = new double[a][b][c];
        for(int i =0; i< a; i++)
		{
			for(int j =0; j< b; j++)
			{
				for(int k =0; k< c; k++)
				{
					Resultado[i][j][k] = 0;
				}
			}
		}
        return Resultado;
    }


    public double[] x_unitario(int dimension)
    {
        double[] x_unitario = new double[dimension];
        for(int i =0; i< dimension; i++) x_unitario[i] = Delta_Kronecker(i, 0);
        return x_unitario;
    }
    
    public double[] y_unitario(int dimension)
    {
        double[] y_unitario = new double[dimension];
        for(int i =0; i< dimension; i++) y_unitario[i] = Delta_Kronecker(i, 1);
        return y_unitario;
    }

    public double[] z_unitario(int dimension)
    {
        double[] z_unitario = new double[dimension];
        for(int i =0; i< dimension; i++) z_unitario[i] = Delta_Kronecker(i, 2);
        return z_unitario;
    }
    
    public double Promedio(double[] vector)
    {
        double suma = 0;
        suma = Suma_Componente_Vector(vector);
        return (suma/vector.length);
    }

	public double[] Promedio_Dos_Vectores(double[] vectorA, double[] vectorB)
	{
		double[] Resultado = new double[vectorA.length];
		if( vectorA.length == vectorB.length)
		{
			for(int i =0; i< vectorA.length; i++)
			{
				Resultado[i] = (vectorA[i] + vectorB[i])*0.5;
			}
		}
		else
		{
			System.out.println("Los vectores no tienen la misma dimension");
		}
		return Resultado;
	}
	
    public double Maximo(int[] vector)
    {
        int maximo = vector[0];
        for(int i=0; i < vector.length; i++) if(maximo < vector[i]) maximo = vector[i];
        return maximo;
    }
    
    public double Maximo(double[] vector)
    {
        double maximo = vector[0];
        for(int i=0; i < vector.length; i++) if(maximo < vector[i]) maximo = vector[i];
        return maximo;
    }

    public double Maximo_Int(List<Integer> lista)
    {
        int maximo = lista.get(0);
        for(int i=0; i < lista.size(); i++) if(maximo < lista.get(i)) maximo = lista.get(i);
        return maximo;
    }
    
    public double Maximo(List<Double> lista)
    {
        double maximo = lista.get(0);
        for(int i=0; i < lista.size(); i++) if(maximo < lista.get(i)) maximo = lista.get(i);
        return maximo;
    }
        
    public int Indice_Maximo(double[] vector)
    {
        double maximo = vector[0];
        int indice = 0;
        for(int i=0; i < vector.length; i++) 
			if(maximo < vector[i]) 
				{
					maximo = vector[i];
					indice = i;
				}
        return indice;
    }

    public double Minimo(int[] vector)
    {
        int minimo = vector[0];
        for(int i=0; i < vector.length; i++) if(minimo > vector[i]) minimo = vector[i];
        return minimo;
    }
        
    public double Minimo(double[] vector)
    {
        double minimo = vector[0];
        for(int i=0; i < vector.length; i++) if(minimo > vector[i]) minimo = vector[i];
        return minimo;
    }

    public double Minimo_Int(List<Integer> lista)
    {
        int minimo = lista.get(0);
        for(int i=0; i < lista.size(); i++) if(minimo < lista.get(i)) minimo = lista.get(i);
        return minimo;
    }
    
    public double Minimo(List<Double> lista)
    {
        double minimo = lista.get(0);
        for(int i=0; i < lista.size(); i++) if(minimo < lista.get(i)) minimo = lista.get(i);
        return minimo;
    }
        
    public int Indice_Minimo(double[] vector)
    {
		
		double minimo = vector[0];
		int indice = 0;
		
		for(int i=0; i < vector.length; i++) 
			if(minimo > vector[i])
			{
				minimo = vector[i];
				indice = i;
			}
        return indice;
    }

    public int Indice_Minimo(List<Double> lista)
    {
		int indice = 0;
		
		if(lista.size() > 0)
		{
			double minimo = lista.get(0);
			
			for(int i=0; i < lista.size(); i++) 
				if(minimo > lista.get(i))
				{
					minimo = lista.get(i);
					indice = i;
				}
		}
        return indice;
    }


	public double[] Dos_Indices_Y_Elementos_Mayores_En_Array_Con_Lista(List<Integer> lista, double[] array)
	{	
		int primer_indice_maximo	= -1; 
		int segundo_indice_maximo 	= -1;
		double primer_valor_maximo 	= Double.NaN;
		double segundo_valor_maximo = Double.NaN;
		
		primer_indice_maximo 	= lista.get(0);
		segundo_indice_maximo 	= -1;
		primer_valor_maximo 	= array[primer_indice_maximo];
		segundo_valor_maximo 	= 0;
			for(int i : lista)
			{
				if(array[i] > primer_valor_maximo)
				{

					segundo_indice_maximo = primer_indice_maximo;
					primer_indice_maximo  = i;

					segundo_valor_maximo = primer_valor_maximo;
					primer_valor_maximo = array[i];
				}
				
				else if(primer_indice_maximo != i && array[i] > segundo_valor_maximo)
				{
					segundo_valor_maximo = array[i];
					segundo_indice_maximo  = i;
				}
				
			}
		
		return new double[]{primer_indice_maximo, segundo_indice_maximo, primer_valor_maximo, segundo_valor_maximo};
	}

	//busca los dos primeros indices mas pequeños en una lista. 
	
    public int[] Indice_Minimo_Dos_Primeros(List<Double> lista)
    {
		int indice = 0;
		int indice_anterior = 0;
		
		//si se buscan dos numeros minimos en una lista, entonces tienen que haber
		//al menos dos elementos
		if(lista.size() > 1)
		{
			//define el primer elemento como el minimo
			double minimo = lista.get(0);
			
			//busca en todos los elementos a ver si hay alguno menor,
			// si hay alguno menor entonces se actualiza el minimo
			// y se guarda el indice en el cual esta ese minimo,
			// el valor anterior del indice se guarda
			for(int i=0; i < lista.size(); i++) 
				if(minimo >= lista.get(i))
				{
					minimo = lista.get(i);
					indice = i;
				}
			if(indice == -1) indice= 0;
			//se reinicia el minimo para buscar el segundo elemento menor			
			minimo = lista.get(0);
			
			for(int i=0; i < lista.size(); i++)
			{ 
				//el segundo menor no puede ser igual al menor
				if(minimo >= lista.get(i) && i != indice)
				{
					minimo = lista.get(i);
					indice_anterior = i;
				}
				
				//el algoritmo tiene un problema cuando justamente el primer
				//elemento es el minimo: ignora todos los condicionales
				//y por lo tanto nunca actualiza el valor por defecto
				//si pasa esto, el valor de "indice" esta bien, el primer
				//elemento es el cero, pero como nunca se entra en los 
				//condicionales, "indice anterior" tambien queda en su
				// valor por defecto, por lo que manualmente hay que pasar
				//al siguiente elemento
				if(indice ==0 && indice_anterior==0)
					minimo = lista.get(1);
			}
			if(indice_anterior == -1) indice_anterior= 0;
		}
		else
			System.out.println("La lista es de una tamanio menor a 2, por lo que no se pueden buscar dos elementos minimos");
			
        return new int[]{indice, indice_anterior};
    }
        
    public double Maximo_2_numeros(double a, double b)
    {
        if(a>=b) return a;
		else return b;
    }
    
    public double Minimo_2_numeros(double a, double b)
    {
        if(a<=b) return a;
		else return b;
    }

	public double Minimo_3_numeros(double a, double b, double c)
    {
		double temporal_1 = 0;
		double temporal_2 = 0;
        if(a<=b) temporal_1 = a;
		if(a>b)  temporal_1 = b;
		if(a<=c) temporal_2 = a;
		if(a>c)  temporal_2 = c;
		if(temporal_1<=temporal_2) 
			return temporal_1;
		else 
			return temporal_2;
    }

	public double Maximo_3_numeros(double a, double b, double c)
    {
		double temporal_1 = 0;
		double temporal_2 = 0;
        if(a>=b) temporal_1 = a;
		if(a<b)  temporal_1 = b;
		if(a>=c) temporal_2 = a;
		if(a<c)  temporal_2 = c;
		if(temporal_1>=temporal_2) 
			return temporal_1;
		else 
			return temporal_2;
    }
	public double Modulo_Vector(double[] vector)
    {
		double suma = 0;
        for(int i=0; i < vector.length; i++) suma += vector[i]*vector[i];
        return Math.sqrt(suma);
    }
    
    public double Angulo_2D(double[] vector)
    {
		return Math.atan2(vector[1], vector[0]);
    }
    
    public double Angulo_2D(double x, double y)
    {
		return Math.atan2(y, x);
    }
    
    
    public double[][] Matriz_En_Vecinos(int[] vecinos, double[][] matriz)
    {
		int k = 0 ;
		double[][] resultado = new double[vecinos.length][matriz[0].length];
		
		for(int i : vecinos)
		{
			resultado[k] = Copiar_Vector(matriz[i]);
			k++;
		}
		
		return resultado;
	}
    
    public double Angulo_2D(double[] vector1, double[] vector2)
    {
		double ang_1 = Angulo_2D(vector1);
		double ang_2 = Angulo_2D(vector2);
		double PI = Param.PI;
		double angulo = ang_1 - ang_2;
		if(ang_1 > 0.5*PI && ang_2 < -0.5*PI)
			angulo  = ang_1 + ang_2;
		if(angulo > PI)
			angulo -= PI;
		if(angulo < -PI)
			angulo += PI;	
		return angulo;
    }
    
    public double[] Vector_Unitario(double[] vector)
    {
		double modulo = Modulo_Vector(vector);
		double UNO_sobre_modulo;
		if(modulo != 0) UNO_sobre_modulo = 1/modulo;
		else   			UNO_sobre_modulo = 0;		
		double[] vector_unitario = new double[vector.length];
		for(int i=0; i < vector.length; i++)  vector_unitario[i] = vector[i]*UNO_sobre_modulo;
		return vector_unitario;
	}

	public double[] Modulo_Vector_Hasta_N(double[][] vector, int N)
    {
		double[] suma = new double[N];
        for(int i=0; i < N; i++)
		{
			for(int j=0; j < vector[0].length; j++)  
			{
				suma[i] += vector[i][j]*vector[i][j];
			}
			suma[i] = Math.sqrt(suma[i]);
		}
        return suma;
    }

    
    public double Delta_Kronecker(double a, double b)
    {
        if(a==b) return 1;
        else return 0;
    }

    public double Delta_Kronecker_Inverso(double a, double b)
    {
        if(a==b) return 0;
        else return 1;
    }

    public double Heavi_Side(double a, double b)
    {
        if(a > b) return 1;
        else return 0;
    }

    public double Varianza(double[] valor_numerico, double valor_real)
    {
        double var = 0;
		int numero_datos = valor_numerico.length;
        for(int i=0; i < numero_datos; i++) var += (valor_numerico[i] - valor_real)*(valor_numerico[i] - valor_real);
        return Math.sqrt(var/numero_datos);
    }
    
    public double Error_Medio_Cuadratico(double[] valor_numerico, double[] valor_real)
    {
        double var = 0;
		int numero_datos = valor_numerico.length;
        for(int i=0; i < numero_datos; i++) var += (valor_numerico[i] - valor_real[i])*(valor_numerico[i] - valor_real[i]);
        return Math.sqrt(var/numero_datos);
    }

    public double Error(double a, double b, double referencia)
    {
        return (a-b)/referencia*100;
    }
    
    public double[] Error_Array(double[] a, double[] b, double[] referencia)
    {
		double[] error = new double[referencia.length];
		for(int i = 0; i<referencia.length; i++)
		{
			if(referencia[i] == 0)
				error[i] = b[i];
			else
				error[i] =  (a[i]-b[i])/referencia[i]*100;
		}
		return error;
    }
    
	public boolean[] Copiar_Vector(boolean[] vector)
	{
		boolean[] Resultado = new boolean[vector.length];
		for(int i =0; i< vector.length; i++)
		{
			Resultado[i] = vector[i];
		}
		return Resultado;
	}
	 
	public int[] Copiar_Vector(int[] vector)
	{
		int[] Resultado = new int[vector.length];
		for(int i =0; i< vector.length; i++)
		{
			Resultado[i] = vector[i];
		}
		return Resultado;
	}

	public double[] Copiar_Vector(double[] vector)
	{
		double[] Resultado = new double[vector.length];
		for(int i =0; i< vector.length; i++)
		{
			Resultado[i] = vector[i];
		}
		return Resultado;
	}

	public List<Integer> Copiar_Lista(List<Integer> lista)
	{
		List<Integer> Resultado = new ArrayList<Integer>();
		for(int i =0; i< lista.size(); i++)
		{
			Resultado.add(lista.get(i));
		}
		return Resultado;
	}

	public List<double[]> Copiar_ListaDD(List<double[]> lista)
	{
		List<double[]> Resultado = new ArrayList<double[]>();
		for(int i =0; i< lista.size(); i++)
		{
			Resultado.add(lista.get(i));
		}
		return Resultado;
	}
	
	// lo que hace es que si se tiene un vector = {1, 2, -8 , 7, valor ... valor}, crea otro vector copiado hasta que empieza a aparecer como entrada la cantidad "valor"
	public double[] Vector_Hasta_Cierto_Valor(double[] vector, double valor)
	{
		double[] Resultado = new double[vector.length];
		for(int i =0; i< vector.length; i++)
		{
			Resultado[i] = vector[i];
		}
		return Resultado;
	}
	
	
	public double[] Copiar_Vector_Hasta_N(double[] vector,int N)
	{
		double[] Resultado = new double[N];
		for(int i =0; i< N; i++)
		{
			Resultado[i] = vector[i];
		}
		return Resultado;
	}
	
	public int[] Copiar_Vector_Hasta_N(int[] vector,int N)
	{
		int[] Resultado = new int[N];
		for(int i =0; i< N; i++)
		{
			Resultado[i] = vector[i];
		}
		return Resultado;
	}

	public int Indice_Distancia_Mas_Corta(int i, int[] lista, double[][] posicion)
	{
		int u=0;
		double[] distancia = new double[lista.length];
		for(int k : lista)
		{
			distancia[u] = Distancia(posicion[k], posicion[i]);
			u++;
		}
		 u = Indice_Minimo(distancia);
		 
		 //if(i == 0) System.out.println(u + " " + lista.length);
		 
		 if(u >= 0)
			return lista[u];
		else
			return -1;
	}
	
	public int Indice_Distancia_Mas_Corta(int i, List<Integer> lista, double[][] posicion)
	{
		int u=0;
		double[] distancia = new double[lista.size()];
		for(int k : lista)
		{
			distancia[u] = Distancia(posicion[k], posicion[i]);
			u++;
		}
		return lista.get(Indice_Minimo(distancia));
	}
	
	public double Distancia(double[] posicion1, double[] posicion2)
	{
		double dx = posicion1[0] - posicion2[0];
		double d2x = dx*dx;
		double dy = posicion1[1] - posicion2[1];
		double d2y = dy*dy;
		return Math.sqrt(d2x + d2y);
	}
	
	public double Pendiente(double[] posicion1, double[] posicion2)
	{
		double dx = posicion1[0] - posicion2[0];
		double dy = posicion1[1] - posicion2[1];
		double pendiente = 0;
		double tolerancia = Math.pow(10,10);
		if(dx == 0 && dy > 0)
			pendiente = tolerancia;
		if(dx == 0 && dy < 0)
			pendiente = - tolerancia;
		if(dx == 0 && dy == 0)
		{
			pendiente = 0;
			System.out.println(" Se esta colocando en el metodo de pendiente dos puntos iguales ");
		}
		if(dx != 0)
			pendiente = dy/dx;	
			
		return pendiente;
	}

    public double[][] Copiar_Matriz_Hasta_N(double[][] matriz, int N)
    {
        double[][] Resultado = new double[N][matriz[0].length];
        for(int i =0; i< N; i++)
        {
            for(int j =0; j< matriz[0].length; j++)
            {
                Resultado[i][j] = matriz[i][j];
            }
        }
        return Resultado;
    }

	public double[][] Copiar_Matriz(double[][] matriz)
	{
		double[][] Resultado = new double[matriz.length][];
		for(int i =0; i< matriz.length; i++)
		{
			Resultado[i]  = new double[matriz[i].length];
			for(int j =0; j< matriz[i].length; j++)
			{
				Resultado[i][j] = matriz[i][j];
			}
		}
		return Resultado;
	}

	
	public int[][] Copiar_Matriz(int[][] matriz)
	{
		int[][] Resultado = new int[matriz.length][];
		for(int i =0; i< matriz.length; i++)
		{
			Resultado[i]  = new int[matriz[i].length];
			for(int j =0; j< matriz[i].length; j++)
			{
				Resultado[i][j] = matriz[i][j];
			}
		}
		return Resultado;
	}

	public boolean[][] Copiar_Matriz(boolean[][] matriz)
	{
		boolean[][] Resultado = new boolean[matriz.length][];
		for(int i =0; i< matriz.length; i++)
		{
			Resultado[i]  = new boolean[matriz[i].length];
			for(int j =0; j< matriz[i].length; j++)
			{
				Resultado[i][j] = matriz[i][j];
			}
		}
		return Resultado;
	}
		public double[][][] Copiar_Tensor_Rango3(double[][][] tensor)
	{
		double[][][] Resultado = new double[tensor.length][tensor[0].length][tensor[0][0].length];
		for(int i =0; i< tensor.length; i++)
		{
			for(int j =0; j< tensor[0].length; j++)
			{
				for(int k =0; k< tensor[0][0].length; k++)
				{
					Resultado[i][j][k] = tensor[i][j][k];
				}
			}
		}
		return Resultado;
	}

    public double[][] Matriz_Rotacion_2D(double angulo_rotacion)
    {
        int dimension = 2;
        double PI_SOBRE_180 = 0.017453292;
        double[][] matriz = new double[dimension][dimension];
        double angulo = PI_SOBRE_180*angulo_rotacion;
        double cos = Math.cos(angulo);
        double sen = Math.sin(angulo);
        matriz[0][0] =   cos;
        matriz[0][1] =   sen;
        matriz[1][0] = - sen;
        matriz[1][1] =   cos;
        return matriz;
    }
    
    public double[][] Derivada_Matriz_Rotacion_2D(double angulo_rotacion)
    {
        int dimension = 2;
        double PI_SOBRE_180 = 0.017453292;
        double[][] matriz = new double[dimension][dimension];
        double angulo = PI_SOBRE_180*angulo_rotacion;
        double cos = Math.cos(angulo);
        double sen = Math.sin(angulo);
        matriz[0][0] = - sen;
        matriz[0][1] =   cos;
        matriz[1][0] = - cos;
        matriz[1][1] = - sen;
        return matriz;
    }

    public double[][] Matriz_Cizalla_2D(double[] cizalla)
    {
        int dimension = 2;
        double[][] matriz = new double[dimension][dimension];
        double PI_SOBRE_180 = 0.017453292;
        double angulo_x = PI_SOBRE_180*cizalla[0];
        double angulo_y = PI_SOBRE_180*cizalla[1];
        matriz[0][0] = 1;
        matriz[0][1] = Math.tan(angulo_x);
        matriz[1][0] = Math.tan(angulo_y);
        matriz[1][1] = 1;
        return matriz;
    }

    public double[][] Matriz_Escado_2D(double[] escala)
    {
        int dimension = 2;
        double[][] matriz = new double[dimension][dimension];
        matriz[0][0] = escala[0];
        matriz[1][1] = escala[1];
        return matriz;
    }

    public double[] Transformacion_Afin_2D(double[] vector, double angulo, double[] vector_traslacion, double[] cizalla, double[] escala)
    {
        // IMPORTANTE:  primero se aplica la escala, luego la cizalla y por ultimo la rotacion
        double[][] matriz_rotacion = Matriz_Rotacion_2D(angulo);
        double[][] matriz_cizalla  = Matriz_Cizalla_2D(cizalla);
        double[][] matriz_escala   = Matriz_Escado_2D(escala);
        double[][] matriz_temporal = Producto_Matrices(matriz_rotacion, matriz_cizalla);
                   matriz_temporal = Producto_Matrices(matriz_temporal, matriz_escala);
        double[] vector_temporal = Producto_Matriz_Vector(matriz_temporal, vector);
        return Suma_Vectores(vector_temporal, vector_traslacion);
    }

    public double[][] Ciclo_Vector_Trans_Afin(double[][] vector_referencia, double angulo, double[] vector_traslacion, double[] cizalla, double[] escala)
    {
        int N = vector_referencia.length;
        int dimension = vector_referencia[0].length;
        double[][] vector = new double[N][dimension];
        for(int i=0; i<N;i++)
        {
            vector[i] = Transformacion_Afin_2D(vector_referencia[i], angulo, vector_traslacion, cizalla,  escala);
        }
        return vector;
    }

    public int N_Total(int[] N)
    {
        int N_Total = 0;
        for(int i =0; i < Param.n_figuras; i++)
        {
            N_Total += N[i];
        }
        return N_Total;
    }

	public double[] Vector_Perpendicular_Horario(double[] vector)
	{
		double[] perpendicular = new double[vector.length];
		perpendicular[0] = vector[1];
		perpendicular[1] = - vector[0];
		
		return perpendicular;
	}

	public double[] Vector_Perpendicular_AntiHorario(double[] vector)
	{
		double[] perpendicular = new double[vector.length];
		perpendicular[0] = - vector[1];
		perpendicular[1] = vector[0];
		
		return perpendicular;
	}
	
	public double[][][] Inv_indices_Matriz(double[][][] matriz)
	{
		double[][][] Resultado = new double[matriz[0].length][matriz[0][0].length][matriz.length];
		for(int i =0; i< matriz[0].length; i++)
		{
			for(int j =0; j< matriz[0][0].length; j++)
			{
				for(int k =0; k< matriz.length; k++)
				{
					Resultado[i][j][k] = matriz[k][i][j];
				}
			}
		}
		return Resultado;
	}


    public double[] Autovalores_Reales_Matriz_2D(double[][] matriz)
    {
		int dimension = matriz.length;
		double[] Resultado = new double[dimension + 1];
		if(dimension !=  matriz[0].length || dimension != 2)
			System.out.println("Este metodo solo calcula los autovalores de matrices cuadradas de dos dimensiones");
		else
		{      
			double traza_2 = 0.5*Traza_R2(matriz);
			double determinante = Determinante_R2_2x2_3x3(matriz);
			double valor_parcial_1 = traza_2*traza_2 - determinante;
			double valor_parcial = 0;
			if(valor_parcial_1 >= 0)
				valor_parcial = Math.sqrt(valor_parcial_1);
			Resultado[0] = traza_2 + valor_parcial;
			Resultado[1] = traza_2 - valor_parcial;

		}
        return Resultado;
    }
    
    public double[] Autovalores_Imaginarios_Matriz_2D(double[][] matriz)
    {
		int dimension = matriz.length;
		double[] Resultado = new double[dimension];
		if(dimension !=  matriz[0].length || dimension != 2)
			System.out.println("Este metodo solo calcula los autovalores de matrices cuadradas de dos dimensiones");
		else
		{      
			double traza_2 = 0.5*Traza_R2(matriz);
			double determinante = Determinante_R2_2x2_3x3(matriz);
			double valor_parcial_1 = traza_2*traza_2 - determinante;
			double valor_parcial = 0;
			if(valor_parcial_1 <= 0)
				valor_parcial = Math.sqrt(-valor_parcial_1);
			Resultado[0] = valor_parcial;
			Resultado[1] = -valor_parcial;

		}
        return Resultado;
    }
    
    public double[] Vector_Componentes_Iguales(int N, double constante)
    {
        double[] Resultado = new double[N];
        for(int i = 0; i < N; i++)
        {
            Resultado[i] = constante;
        }
        return Resultado;
    }
    
    public int[] Vector_Componentes_Iguales(int N, int constante)
    {
        int[] Resultado = new int[N];
        for(int i = 0; i < N; i++)
        {
            Resultado[i] = constante;
        }
        return Resultado;
    }
    
    public boolean[] Vector_Componentes_Iguales(int N, boolean constante)
    {
        boolean[] Resultado = new boolean[N];
        for(int i = 0; i < N; i++)
        {
            Resultado[i] = constante;
        }
        return Resultado;
    }
    
    public double[][] Matriz_Componentes_Iguales(int N, int M,  double constante)
    {
        double[][] Resultado = new double[N][M];
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < M; j++)
			{
				Resultado[i][j] = constante;
			}
        }
        return Resultado;
    }

    public int[][] Matriz_Componentes_Iguales(int N, int M,  int constante)
    {
        int[][] Resultado = new int[N][M];
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < M; j++)
			{
				Resultado[i][j] = constante;
			}
        }
        return Resultado;
    }
        
    public double[][] Matriz_Vectores_Iguales(int N, double[] constante)
    {
        double[][] Resultado = new double[N][constante.length];
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < constante.length; j++)
			{
				Resultado[i][j] = constante[j];
			}
        }
        return Resultado;
    }
    
    public double[] Vector_Centro_de_Masas(double[][] vector)
    {
		int N 		  = vector.length;
		int dimension = vector[0].length;     
		double[] vector_centro_de_masa = new double[dimension];
		
		for(int i = 0; i< N; i++)
        {
			vector_centro_de_masa[0] += vector[i][0];
			vector_centro_de_masa[1] += vector[i][1];							
        }
        
		vector_centro_de_masa[0] = 1/(double)N*vector_centro_de_masa[0];	
		vector_centro_de_masa[1] = 1/(double)N*vector_centro_de_masa[1];

		
		return vector_centro_de_masa;
	}

public int[] Siguiente_Elemento_Sentido_Horario(int[][] vecinos_borde, double[][] posicion)
	{
		int 	largo 			= vecinos_borde[0].length;
		int 	k 				= 0;
		int[] 	sentido_horario = new int[largo];
				
		double 	momento;
		double 	dx;
		double 	dy;

		
		for(int i : vecinos_borde[0])
		{		
			//System.out.println( i + " " + largo + " " + vecinos_borde[1][k]);
			dx = posicion[i][0] - posicion[vecinos_borde[1][k]][0];
			dy = posicion[i][1] - posicion[vecinos_borde[1][k]][1];
						
			momento = Producto_Cruz(new double[] {dx,dy}, posicion[i])[0];
			
			if(momento > 0)
				sentido_horario[k] = vecinos_borde[1][k];
			else
				sentido_horario[k] = vecinos_borde[2][k];
			
			//System.out.println(i+ " " + sentido_horario[k]	);
			k++;
		}
		return sentido_horario;
	}
	
	public int[] Siguiente_Elemento_Sentido_AntiHorario(int[][] vecinos_borde, double[][] posicion)
	{
		int 	largo 			= vecinos_borde[0].length;
		int 	k 				= 0;
		int[] 	sentido_horario = new int[largo];
				
		double 	momento;
		double 	dx;
		double 	dy;

		
		for(int i : vecinos_borde[0])
		{		
			
			dx = posicion[i][0] - posicion[vecinos_borde[1][k]][0];
			dy = posicion[i][1] - posicion[vecinos_borde[1][k]][1];
						
			momento = Producto_Cruz(new double[] {dx,dy}, posicion[i])[0];
			
			if(momento < 0)
				sentido_horario[k] = vecinos_borde[1][k];
			else
				sentido_horario[k] = vecinos_borde[2][k];
			
			//System.out.println(i+ " " + sentido_horario[k]	);
			k++;
		}
		return sentido_horario;
	}
	
	
	public double[] Reordenar_Vector(double[] vector, int[] lista)
	{
		int k = 0;
		double[] Resultado = new double[vector.length];
		for(int i: lista)
		{
			Resultado[k] =  vector[i];
			k++;
		}	
		return Resultado;
	}

	public int[] Reordenar_Vector(int[] vector, int[] lista)
	{
		int k = 0;
		int[] Resultado = new int[vector.length];
		for(int i: lista)
		{
			Resultado[k] =  vector[i];
			k++;
		}	
		return Resultado;
	}
	
	public int[] Recortar_Vector(int[] vector, int indice_inicial, int indice_final)
	{
		int k = 0;
		int[] Resultado = new int[indice_final - indice_inicial];
		for(int i = indice_inicial; i < indice_final; i++)
		{
			Resultado[k] =  vector[i];
			//System.out.println(i + " " + Resultado.length );
			k++;
		}	
		return Resultado;
	}
	
	public double[] Recortar_Vector(double[] vector, int indice_inicial, int indice_final)
	{
		int k = 0;
		double[] Resultado = new double[indice_final - indice_inicial];
		for(int i = indice_inicial; i < indice_final; i++)
		{
			Resultado[k] =  vector[i];
			//System.out.println(i + " " + Resultado.length );
			k++;
		}	
		return Resultado;
	}
	
    public  void Escribir_en_Tiempo(double[] variables, String[] nombres)
	{
		String directorio_serie_temporal = Param.serie_temporal;
		File serie_temporal 		 	 = new File(directorio_serie_temporal);

		try
		{
			FileWriter archivo = new FileWriter(serie_temporal,true);
			PrintWriter escribir = new PrintWriter(archivo);
            for(int entradas = 0; entradas < variables.length; entradas++)
            {
                 escribir.write(" " + variables[entradas]);
            }
			escribir.write(System.getProperty("line.separator"));
			escribir.close();
		}  
		catch (IOException a) { System.err.println("Error al escribir");}

    }
    
    public  void Escribir_en_Espacio(double[][] variables, String[] nombres, int contador_tiempo)
	{
		int numero_variables  = variables.length;
		int numero_particulas = variables[0].length;
		String directorio_serie_espacial = Param.serie_espacial + "_" + String.valueOf(contador_tiempo) + ".txt";
		File serie_espacial 			 = new File(directorio_serie_espacial);

		try
		{
			FileWriter archivo = new FileWriter(serie_espacial,true);
			PrintWriter escribir = new PrintWriter(archivo);
			for(int i = 0; i < numero_particulas; i++)
            {
				for(int entradas = 0; entradas < numero_variables; entradas++)
				{
					 escribir.write(" " + variables[entradas][i]);
				}
			escribir.write(System.getProperty("line.separator"));	
			}	
			escribir.close();
		}  
		catch (IOException a) { System.err.println("Error al escribir");}

    }
    
    
    public void Limpiar_Archivos(File carpeta, String extension)
	{
		File fList[] = carpeta.listFiles();
		for (File f : fList) 
		{
			if (f.getName().endsWith("." + extension)) 
			{
				f.delete(); 
			}
		}
	}
	
	public boolean Pertenece_A_La_Lista(int i, int[] lista)
	{
		boolean pertenece = false;
		int N = lista.length;
		
		for(int j : lista)
		{
			if(i == j)
			pertenece = true;
		}
		
		return pertenece;
	}	
	
	

	public void Imprimir_Comparar_Vector(int[] vectorA, int[] vectorB)
	{
		System.out.println(" ");
		if(vectorA.length != vectorB.length)
		{
			System.out.println(" Los vectores que se estan comparando en el metodo Imprimir_Comparar_Vector (int) no tienen la misma longitud");
			Limpiar_Archivos(new File(Param.directorio_base), "class");
			System.exit(1);
		}
		else
		{
			for(int i = 0; i < vectorA.length; i++)
			{
				System.out.println(i + " " + vectorA[i] + " " + vectorB[i]);
			}
			System.out.println(" ");
		}
	}

	public void Imprimir_Comparar_Vector(int[] vectorA, double[] vectorB)
	{
		System.out.println(" ");
		if(vectorA.length != vectorB.length)
		{
			System.out.println(" Los vectores que se estan comparando en el metodo Imprimir_Comparar_Vector (int, double) no tienen la misma longitud");
			Limpiar_Archivos(new File(Param.directorio_base), "class");
			System.exit(1);
		}
		else
		{
			for(int i = 0; i < vectorA.length; i++)
			{
				System.out.println(i + " " + vectorA[i] + " " + vectorB[i]);
			}
			System.out.println(" ");
		}
	}
	
	public void Imprimir_Comparar_Vector(double[] vectorA, double[] vectorB)
	{
		System.out.println(" ");
		if(vectorA.length != vectorB.length)
		{
			System.out.println(" Los vectores que se estan comparando en el metodo Imprimir_Comparar_Vector (double) no tienen la misma longitud");
			Limpiar_Archivos(new File(Param.directorio_base), "class");
			System.exit(1);
		}
		else
		{
			for(int i = 0; i < vectorA.length; i++)
			{
				System.out.println(i + " " + vectorA[i] + " " + vectorB[i]);
			}
			System.out.println(" ");
		}
	}
	
	public void Imprimir_Comparar_Vector(boolean[] vectorA, boolean[] vectorB)
	{
		System.out.println(" ");
		if(vectorA.length != vectorB.length)
		{
			System.out.println(" Los vectores que se estan comparando en el metodo Imprimir_Comparar_Vector (boolean) no tienen la misma longitud");
			Limpiar_Archivos(new File(Param.directorio_base), "class");
			System.exit(1);
		}
		else
		{
			for(int i = 0; i < vectorA.length; i++)
			{
				System.out.println(i + " " + vectorA[i] + " " + vectorB[i]);
			}
			System.out.println(" ");
		}
	}
	
	public void Imprimir_Lista(List<Integer> Lista)
	{
		System.out.println(Arrays.toString(Lista.toArray()));
	}


	public void Imprimir_ListaD(List<Double> Lista)
	{
		System.out.println(Arrays.toString(Lista.toArray()));
	}
	
	public void Imprimir_ListaDD(List<double[]> Lista)
	{
		for(int i = 0; i < Lista.size(); i++)
		System.out.println("[" + Lista.get(i)[0] + "," +  Lista.get(i)[1] + "]");
	}


	public void Imprimir_Vector_Transpuesto(int[] vector)
	{
		System.out.println(" ");
		for(int i = 0; i < vector.length; i++)
		{
			System.out.print(i + " ");
		}
		System.out.println(" ");
		for(int i = 0; i < vector.length - 1; i++)
		{
			System.out.print(vector[i] + ", ");
		}
		System.out.print(vector[vector.length - 1]);	
		System.out.println(" ");	
	}
	
	public void Imprimir_Vector_Transpuesto(double[] vector)
	{
		System.out.println(" ");
		for(int i = 0; i < vector.length; i++)
		{
			System.out.print(i + " ");
		}
		System.out.println(" ");
		for(int i = 0; i < vector.length; i++)
		{
			System.out.print(vector[i] + ", ");
		}
		System.out.print(vector[vector.length - 1]);
		System.out.println(" ");		
	}
	
	public void Imprimir_Vector_Transpuesto(boolean[] vector)
	{
		System.out.println(" ");
		for(int i = 0; i < vector.length; i++)
		{
			System.out.print(i + " ");
		}
		System.out.println(" ");
		for(int i = 0; i < vector.length; i++)
		{
			System.out.print(vector[i] + ", ");
		}
		System.out.print(vector[vector.length - 1]);
		System.out.println(" ");		
	}
		
	
	
	public void Imprimir_Vector(int[] vector)
	{
		System.out.println(" ");
		for(int i = 0; i < vector.length; i++)
		{
			System.out.println(i + " " + vector[i]);
		}	
		System.out.println(" ");	
	}
	
	public void Imprimir_Vector(int[] vector, int[] lista)
	{
		System.out.println(" ");
		int k = 0;
		for(int i: lista)
		{
			System.out.println(k + " " + vector[i]);
			k++;
		}	
		System.out.println(" ");	
	}
	
	public void Imprimir_Vector(double[] vector, int[] lista)
	{
		System.out.println(" ");
		int k = 0;
		for(int i: lista)
		{
			System.out.println(k + " " + vector[i]);
			k++;
		}	
		System.out.println(" ");	
	}
	
	public void Imprimir_Vector(boolean[] vector, int[] lista)
	{
		System.out.println(" ");
		int k = 0;
		for(int i: lista)
		{
			System.out.println(k + " " + vector[i]);
			k++;
		}	
		System.out.println(" ");	
	}
	
	public void Imprimir_Vector(double[] vector)
	{
		System.out.println(" ");
		for(int i = 0; i < vector.length; i++)
		{
			System.out.println(i + " " + vector[i]);
		}
		System.out.println(" ");		
	}
	
	public void Imprimir_Vector(boolean[] vector)
	{
		System.out.println(" ");
		for(int i = 0; i < vector.length; i++)
		{
			System.out.println(i + " " + vector[i]);
		}
		System.out.println(" ");		
	}
	
	public void Imprimir_Matriz(int[][] matriz)
	{
		System.out.println(" ");
		for(int i = 0; i < matriz.length; i++)
		{
			System.out.print(i + " ");
			for(int j = 0; j < matriz[i].length - 1; j++)
				System.out.print(matriz[i][j] + ", ");
			if(matriz[i].length > 1)System.out.print(matriz[i][matriz[i].length - 1] + " ");	
			System.out.println(" ");
		}
		System.out.println(" ");
	}		
	
	public void Imprimir_Matriz(boolean[][] matriz)
	{
		System.out.println(" ");
		for(int i = 0; i < matriz.length; i++)
		{
			System.out.print(i + " ");
			for(int j = 0; j < matriz[i].length - 1; j++)
				System.out.print(matriz[i][j] + ", ");
			if(matriz[i].length > 1)System.out.print(matriz[i][matriz[i].length - 1] + " ");	
			System.out.println(" ");
		}
		System.out.println(" ");
	}	
	
	public void Imprimir_Matriz(double[][] matriz)
	{
		System.out.println(" ");
		for(int i = 0; i < matriz.length; i++)
		{
			System.out.print(i + " ");
			for(int j = 0; j < matriz[i].length - 1; j++)
				System.out.print(matriz[i][j] + ", ");
			if(matriz[i].length > 1)System.out.print(matriz[i][matriz[i].length - 1] + " ");	
			System.out.println(" ");
		}
		System.out.println(" ");
	}
	
		public void Imprimir_Matriz_Transpuesta(int[][] matriz)
	{
		System.out.println(" ");
		for(int i = 0; i < matriz[0].length; i++)
		{
			for(int j = 0; j < matriz.length; j++)
				System.out.print(matriz[j][i] + " ");
			System.out.println(" ");
		}
		System.out.println(" ");
	}		
	
	public void Imprimir_Matriz_Transpuesta(boolean[][] matriz)
	{
		System.out.println(" ");
		for(int i = 0; i < matriz[0].length; i++)
		{
			for(int j = 0; j < matriz.length; j++)
				System.out.print(matriz[j][i] + " ");
			System.out.println(" ");
		}
		System.out.println(" ");
	}	
	
	public void Imprimir_Matriz_Transpuesta(double[][] matriz)
	{
		System.out.println(" ");
		for(int i = 0; i < matriz[0].length; i++)
		{
			for(int j = 0; j < matriz.length; j++)
				System.out.print(matriz[j][i] + " ");
			System.out.println(" ");
		}
		System.out.println(" ");
	}
	
	public int[] Lista_Primer_Indice_Sparce_Matrix(int[][] sparce_matrix)
	{
		int N = sparce_matrix.length;
		int[] lista = new int[N];
		
		for(int i = 0 ; i < N; i++)
		{
			lista[i] = sparce_matrix[i][0];
		}
		
		return lista;
	}
	
	public double[][] Juntar_Vectores(double[] vector_x, double[] vector_y)
	{
		double[][] resultado = new double[vector_x.length][2];
		for(int i = 0 ; i < vector_x.length; i++)
		{
			resultado[i][0] = vector_x[i];
			resultado[i][1] = vector_y[i];
		}
		return resultado;
	}
	

	/*
	//Este metodo lo que hace es buscar en una lista un elemento y sacarlo de esa lista
	public int[] Sacar_Elemento_Lista(int elemento, int[] lista)
	{
		int k = 0;
		int[] nueva_lista = new int[lista.length - 1];
		for(int i : lista)
        {
			if(i != elemento)
			{
				
				nueva_lista[k] = i;
				k++;
			}
			if(i == 71)System.out.println(elemento);
        }
        return nueva_lista;
	}
	*/
	
	public static int[] Sacar_Elemento_Lista(int elemento, int[] lista)
	{
		List<Integer> resultado  	= new ArrayList<Integer>();

		for(int i : lista)
			if(i != elemento)
				resultado.add(i);

		return resultado.stream().filter(ñ -> ñ != null).mapToInt(ñ -> ñ).toArray();	
	}
	
	
	public static double[] Sacar_Elemento_Lista(double elemento, double[] lista)
	{
		List<Double> resultado  	= new ArrayList<Double>();

		for(int i = 0; i < lista.length; i++)
			if(lista[i] != elemento)
				resultado.add(lista[i]);

		return resultado.stream().filter(ñ -> ñ != null).mapToDouble(ñ -> ñ).toArray();	
	}

	//este metodo busca el conjuntos de elementos viejos en la listay lo reemplaza por el elemento nuevo
	public int[] Reemplazar_Elementos_Lista(int[] elementos_viejos, int elemento_nuevo, int[] lista)
	{
		int k = 0;
		int[] nueva_lista = new int[lista.length];
		for(int i : lista)
        {
			if(Indice_De_Elemento(i, elementos_viejos) == -1)
			{
				nueva_lista[k] = i;
				
			}
			else
				nueva_lista[k] = elemento_nuevo;
		k++;		
        }
        return nueva_lista;
	}
	
	//este metodo busca el elemento_viejo en la lista y lo cambia por el elemento_nuevo
	public int[] Reemplazar_Elemento_Lista(int elemento_viejo, int elemento_nuevo, int[] lista)
	{
		int k = 0;
		int[] nueva_lista = new int[lista.length];
		for(int i : lista)
        {
			if(i != elemento_viejo)
			{
				nueva_lista[k] = i;
				
			}
			else
				nueva_lista[k] = elemento_nuevo;
		k++;		
        }
        return nueva_lista;
	}

	//este metodo busca el elemento_viejo en la lista y lo cambia por el elemento_nuevo
	public double[] Reemplazar_Elemento_Lista(double elemento_viejo, double elemento_nuevo, double[] lista)
	{
		int k = 0;
		double[] nueva_lista = new double[lista.length];
		
		for(int i = 0; i < lista.length; i++)
        {
			if(lista[i] != elemento_viejo)
			{
				nueva_lista[k] = lista[i];
				
			}
			else
				nueva_lista[k] = elemento_nuevo;
		k++;		
        }
        return nueva_lista;
	}
	
	public int Siguiente_Elemento_Diferente_De(int elemento_no_deseado, int[] lista)
	{
		int k = elemento_no_deseado;
		for(int i : lista)
		{
			if(i != elemento_no_deseado)
			{
				k = i;
				break;
			}		
		}
		return k;	
	}
	
	public int Indice_De_Elemento(double elemento, double[] lista)
	{
		int indice = -1;
		for(int i = 0; i < lista.length; i++)
		{
			if(elemento == lista[i])
				indice = i;
		}
		return indice;
	}
	
	public List<Integer> Lista_Indice_De_Elemento(double elemento, double[] lista)
	{
		List<Integer> indice = new ArrayList<Integer>();
		for(int i = 0; i < lista.length; i++)
		{
			if(elemento == lista[i])
				indice.add(i);
		}
		return indice;
	}
	
	public int Indice_De_Elemento(double elemento, double[][] lista)
	{
		int indice = -1;
		for(int i = 0; i < lista.length; i++)
		{
			if(elemento == lista[i][0])
				indice = i;
		}
		return indice;
	}
	
	public int Indice_De_Elemento(int elemento, int[][] lista)
	{
		int indice = -1;
		for(int i = 0; i < lista.length; i++)
		{
			if(elemento == lista[i][0])
				indice = i;
		}
		return indice;
	}
	
	public int Indice_De_Elemento(int elemento, int[] lista)
	{
		int indice = -1;
		int k = 0;
		if(lista.length > 0)
		{
			for(int j : lista)
			{
				if(elemento == j)
					indice = k;
				k++;
			}
		}
		return indice;
	}
	
	public double[] Ordenar_Vector(double[] vector)
	{
		int largo = vector.length;
		for (int i = 0; i < largo; i++)
        {
            int index = i;
            for (int j = i + 1; j < largo; j++)
                if (vector[j] < vector[index]) 
                    index = j;
      
            double smallerNumber = vector[index];  
            vector[index] = vector[i];
            vector[i] = smallerNumber;
        }
        return vector;
	}
}

	
class  Ordenar_Vector_Menor_Mayor
{
	Utilidades Util = new Utilidades();
	int[] indice;
	int[] array_int;
	double[] array_double;
	
	public Ordenar_Vector_Menor_Mayor(double[] vector)
	{
		//int k = 0;
		int largo = vector.length;
		int indice_actual = 0;
		int indice_anterior = -1;
		array_double = Util.Copiar_Vector(vector);
		
		indice = new int[largo];
		List<Integer> resultado  	= new ArrayList<Integer>();
		
		for (int i = 0; i < largo ; i++)
        {
            int index = i;
            //System.out.println(i + " " + k + " " + index + " " + largo);
            for (int j = index + 1; j < vector.length ; j++)
                if (array_double[j] < array_double[index])
                    index = j;
			
            double smallerNumber = array_double[index]; 
        
			array_double[index] = array_double[i];
			array_double[i] 	= smallerNumber;
			indice[i] 			= Util.Indice_De_Elemento(smallerNumber,vector);
			//Util.Imprimir_Lista( Util.Lista_Indice_De_Elemento(smallerNumber, vector));
			//System.out.println(i + " " + index);
			

          //  System.out.println(smallerNumber + " " + i + " " + index);
        }
        
		//indice = resultado.stream().filter(ñ ->ñ != null).mapToInt(ñ -> ñ).toArray();

	}
	public Ordenar_Vector_Menor_Mayor(int[] vector)
	{
		int k = 0;
		int largo = vector.length;
		int indice_actual = 0;
		int indice_anterior = -1;
		array_int = Util.Copiar_Vector(vector);
		
		indice = new int[largo];

		List<Integer> resultado  	= new ArrayList<Integer>();
		
		for (int i = 0; i < largo ; i++)
        {
            int index = k;
            for (int j = i + 1; j < largo ; j++)
                if (array_int[j] < array_int[index])
                { 
                    index = j;
                    
                 }   
      
			
            int smallerNumber = array_int[index]; 
            indice_anterior = indice_actual;
            indice_actual = Util.Indice_De_Elemento(smallerNumber,vector);
			
            if(i > 0 && indice_actual == indice_anterior)
            {
				i++;
				largo++;
			}
            else
            {
				
				resultado.add(indice_actual);
				array_int[index] = array_int[k];
				array_int[k] = smallerNumber;
				k++;
			}
            // System.out.println(i + " " + index + " " + smallerNumber + " " + resultado[0][i] + " " + resultado[0][index] + " " + Indice_de_Elemento(smallerNumber,vector[0]));

            

          //  System.out.println(smallerNumber + " " + i + " " + index);
        }
		indice = resultado.stream().filter(ñ ->ñ != null).mapToInt(ñ -> ñ).toArray();
	}
}

class  Ordenar_Vector
{
	Utilidades Util = new Utilidades();
	int[] indices;
	int[] array_int;
	double[] array_double;
	

	public Ordenar_Vector(double[] vector)
	{
		int index 		= 0;
		int largo 		= vector.length;
		indices 		= new int[largo];
		array_double 	= new double[largo];
				
		HashMap<Integer, Double> map = new HashMap<Integer, Double>();
		
		for (int i = 0; i < largo ; i++)
			map.put(i, vector[i]);
		
		Map<Integer, Double> sortedMap = map.entrySet().stream().sorted(Entry.comparingByValue()).collect(Collectors.toMap(Entry::getKey, Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new));
              

		for (Map.Entry<Integer, Double> mapEntry : sortedMap.entrySet()) 
		{
			indices[index] 		= mapEntry.getKey();
			array_double[index] = mapEntry.getValue();
			index++;
		}  
	}
	
	public Ordenar_Vector(int[] vector)
	{
		int index 		= 0;
		int largo 		= vector.length;
		indices 		= new int[largo];
		array_int 	= new int[largo];
				
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		
		for (int i = 0; i < largo ; i++)
			map.put(i, vector[i]);
		
		Map<Integer, Integer> sortedMap = map.entrySet().stream().sorted(Entry.comparingByValue()).collect(Collectors.toMap(Entry::getKey, Entry::getValue, (e1, e2) -> e1, LinkedHashMap::new));
              

		for (Map.Entry<Integer, Integer> mapEntry : sortedMap.entrySet()) 
		{
			indices[index] 		= mapEntry.getKey();
			array_double[index] = mapEntry.getValue();
			index++;
		}  
	}
}


	
