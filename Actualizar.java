public class Actualizar
{
    public double Tensor_R0(double escalar_anterior, double tasa_de_cambio_escalar, double tiempo )
    {
        return escalar_anterior + tasa_de_cambio_escalar*tiempo;
    }

    public double[] Tensor_R1(double[] vector_anterior, double[] tasa_de_cambio_vector, double tiempo )
    {
        double[] Resultado = new double[vector_anterior.length];
        for(int i = 0; i< vector_anterior.length; i++)
        {
            Resultado[i] = vector_anterior[i] + tasa_de_cambio_vector[i]*tiempo;
        }
        return Resultado;
    }

    public double[][] Tensor_R2(double[][] matriz_anterior, double[][] tasa_de_cambio_matriz, double tiempo )
    {
        double[][] Resultado = new double[matriz_anterior.length][matriz_anterior[0].length];
        for(int i = 0; i< matriz_anterior.length; i++)
        {
            for(int j = 0; j< matriz_anterior[0].length; j++)
            {
                Resultado[i][j] = matriz_anterior[i][j] + tasa_de_cambio_matriz[i][j]*tiempo;
            }
        }
        return Resultado;
    }

    public double Tensor_R0_Final(double escalar, double escalar_anterior )
    {
        return 2.0*escalar - escalar_anterior;
    }

    public double[] Tensor_R1_Final(double[] vector, double[] vector_anterior )
    {
        double[] Resultado = new double[vector_anterior.length];
        for(int i = 0; i< vector_anterior.length; i++)
        {
            Resultado[i] = 2.0*vector[i] - vector_anterior[i];
        }
        return Resultado;
    }

    public double[][] Tensor_R2_Final(double[][] matriz, double[][] matriz_anterior)
    {
        double[][] Resultado = new double[matriz_anterior.length][matriz_anterior[0].length];
        for(int i = 0; i< matriz_anterior.length; i++)
        {
            for(int j = 0; j< matriz_anterior[0].length; j++)
            {
                Resultado[i][j] = 2.0*matriz[i][j]  - matriz_anterior[i][j];
            }
        }

        return Resultado;
    }

    public double Tensor_R0_Final_PC(double escalar_anterior, double tasa_de_cambio_escalar, double tasa_de_cambio_escalar_anterior , double tiempo)
    {
        return escalar_anterior + 0.5*tiempo*(tasa_de_cambio_escalar + tasa_de_cambio_escalar_anterior);
    }
    
    public double[][] Tensor_R2_Final_PC(double[][] matriz_anterior, double[][] tasa_de_cambio_matriz, double[][] tasa_de_cambio_matriz_anterior , double tiempo)
    {
        double[][] Resultado = new double[matriz_anterior.length][matriz_anterior[0].length];
        for(int i = 0; i< matriz_anterior.length; i++)
        {
            for(int j = 0; j< matriz_anterior[0].length; j++)
            {
                Resultado[i][j] = matriz_anterior[i][j] + 0.5*tiempo*(tasa_de_cambio_matriz[i][j] + tasa_de_cambio_matriz_anterior[i][j]);
            }
        }

        return Resultado;
    }

    public double[] Posicion_Gray(double[] vector_anterior, double[] velocidad, double[] aceleracion, double tiempo )
    {
        double[] vector = new double[vector_anterior.length];
        for(int i = 0; i< vector_anterior.length; i++)
        {
            vector[i] = vector_anterior[i] + velocidad[i]*tiempo + 0.5*aceleracion[i]*tiempo*tiempo;
        }
        return vector;
    }
    


    public double Tensor_R0_Final_Gray(double escalar_predicho, double tasa_de_cambio_escalar_predicho, double tasa_de_cambio_escalar_anterior , double tiempo)
    {
        return escalar_predicho + 0.5*tiempo*(tasa_de_cambio_escalar_predicho - tasa_de_cambio_escalar_anterior);
    }

    public double[] Tensor_R1_Final_Gray(double[] vector_predicho, double[] tasa_de_cambio_vector_predicho, double[] tasa_de_cambio_vector_anterior , double tiempo)
    {
        double[] Resultado = new double[vector_predicho.length];
        for(int i = 0; i< vector_predicho.length; i++)
        {
            Resultado[i] =  vector_predicho[i] + 0.5*tiempo*(tasa_de_cambio_vector_predicho[i] - tasa_de_cambio_vector_anterior[i]);
        }
        return Resultado;
    }

    public double[][] Tensor_R2_Final_Gray(double[][] matriz_predicha, double[][] tasa_de_cambio_matriz_predicho, double[][] tasa_de_cambio_matriz_anterior , double tiempo)
    {
        double[][] Resultado = new double[matriz_predicha.length][matriz_predicha[0].length];
        for(int i = 0; i< matriz_predicha.length; i++)
        {
            for(int j = 0; j< matriz_predicha[0].length; j++)
            {
                Resultado[i][j] = matriz_predicha[i][j] + 0.5*tiempo*(tasa_de_cambio_matriz_predicho[i][j] - tasa_de_cambio_matriz_anterior[i][j]);
            }
        }

        return Resultado;
    }

    public double[] Posicion_Mod(double[] vector_anterior, double[] velocidad_anterior, double[] velocidad, double tiempo )
    {
        double[] vector = new double[vector_anterior.length];
        for(int i = 0; i< vector_anterior.length; i++)
        {
            vector[i] = vector_anterior[i] + 0.5*(velocidad_anterior[i] + velocidad[i])*tiempo;
        }
        return vector;
    }
}

