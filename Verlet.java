import java.util.List;

public class Verlet
{

    //objetos de los metodos
	Utilidades Util  	= new Utilidades();
	Parametros Param 	= new Parametros();
    Actualizar Actual 	= new Actualizar();
	Condiciones_Borde	Condicion_Borde = new Condiciones_Borde();


    
	//Variables que se actualizan en cada paso de tiempo
	double[][] posicion;
	double[][] velocidad_intermedia;
	double[]   posicion_centro_de_masa;
	double[][] posicion_particulas_centro_de_masa;
	double[][] elongacion;
	


    public Verlet
    (
		double 	     delta_tiempo_anterior,
		List<Integer>	particulas_en_borde,
		int[]	 	 borde_de_particula, 
		double[][]   posicion_anterior, 
		double[][]   posicion_inicial, 
		double[]	 posicion_centro_de_masa_inicial,
		double[][]   velocidad_anterior, 
		double[][] 	 aceleracion_anterior
	)
    {	
		
		//variables auxiliares
		boolean		esta_particula_en_borde = false;
		int 		N 		  				= posicion_inicial.length;
        int 		dimension 				= posicion_inicial[0].length;
        double 		medio_delta_tiempo 		= 0.5*delta_tiempo_anterior;
        double 		angulo_centro_de_masas 	= 0;
		double[]  	elongacion_sin_borde	= new double[dimension];
        double[] 	angulo_cada_particula 	= new double[N];
        double[] 	angulo_cada_particula2 	= new double[N];
        double[] 	dif_centro_masas 		= new double[dimension];
        double[][]	matriz_rotacion			= new double[dimension][dimension];
        double[][] 	posicion_relajacion 	= new double[N][dimension];
        
		//Variables que se actualizan en cada paso de tiempo
        posicion 			 				= new double[N][dimension];
        velocidad_intermedia 				= new double[N][dimension];
		elongacion 						 	= new double[N][dimension];
		posicion_particulas_centro_de_masa 	= new double[N][dimension];
	
        for(int i = 0; i< N; i++)
        {
            velocidad_intermedia[i]  = Actual.Tensor_R1(velocidad_anterior[i], aceleracion_anterior[i], medio_delta_tiempo);
            posicion[i]    			 = Actual.Tensor_R1(posicion_anterior[i],  velocidad_intermedia[i], delta_tiempo_anterior);
		}
        
		posicion_centro_de_masa = Util.Vector_Centro_de_Masas(posicion);

          for(int i = 0; i< N; i++)
        {
			esta_particula_en_borde 				= borde_de_particula[i] != -2;
            posicion_particulas_centro_de_masa[i] 	= Util.Resta_Vectores(posicion[i], posicion_centro_de_masa);
            elongacion[i]     			  			= Util.Resta_Vectores(posicion[i], posicion_inicial[i]);
			elongacion[i] 						  	= Condicion_Borde.Condiciones_Desplazamiento(i, esta_particula_en_borde, borde_de_particula[i], elongacion[i]);
			posicion[i]  						  	= Util.Suma_Vectores(elongacion[i], posicion_inicial[i]);
        }
        //System.out.println(posicion_inicial[0][0] + " despues " + posicion[0][0]);
    }	

}
