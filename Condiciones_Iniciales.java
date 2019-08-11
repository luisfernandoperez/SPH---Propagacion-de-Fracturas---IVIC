//En este metodo se inicializa aquellas variables del material EN EQUILIBRIO TERMODINAMICO, es decir, en reposon, sin elongaciones ni estreses internos o externos, las variables de SPH son inicializadas aca, 
import java.util.*;
import java.util.stream.*;

public class Condiciones_Iniciales
{
    Utilidades 		Util 		= new Utilidades();
    Otras_Variables Otras_Var 	= new Otras_Variables();
    Parametros 		Param 		= new Parametros();
    Tasa_Variables 	Tasa_Var   	= new Tasa_Variables();
    Deteccion_Bordes Bordes		= new Deteccion_Bordes();
	Actualizar 		Actual 		= new Actualizar();
	Kernel_y_Derivadas Adap_y_Consis = new Kernel_y_Derivadas();
	Particulas_Frontera Frontera = new Particulas_Frontera();
	Condiciones_Borde	Condicion_Borde = new Condiciones_Borde();	
	
	//Variables basicas
	int    dimension = Param.dimension;
	double tiempo;
	double masas;
	double RAIZ_DOS = Param.RAIZ_DOS;


	//Variables propias de SPH
	int		 	 numero_vecinos;
	double[]	 h;
	int[][]		 vecinos;
	int[][] 	 vecinos_potenciales; 	
	int[][] 	 celda_de_particula;
	double[][]	 kernel_sin_consistencia;
	double[][]	 kernel;
    double[][][] grad_kernel_sin_consistencia;
    double[][][] grad_kernel;

	//Variables que se actualizan en cada paso de tiempo
	double[]  	 energia_elastica ;
	double[]	 energia_cinetica;
	double[][]	 posicion;
	double[][] 	 posicion_sin_deformar;
	double[][]	 velocidad;
	double[] 	 densidad;	
	
	//Tasa de cambio de las variables de arriba
	double[] 	 tasa_energia_elastica ;
	double[][]	 aceleracion;
	
    //Variables relacionadas con el centro de masa
	double   	angulo_centro_de_masas;
	double[][] 	posicion_particulas_centro_de_masa;  
	double[] 	posicion_centro_de_masa;
	double[] 	velocidad_centro_de_masa;
	double[][] 	velocidad_particulas_centro_de_masa;
		
	// Variables relacionadas con el calculo de la ecuacion de movimiento
	double[][] 	 elongacion;
    double[][][] tensor_deformacion;
	double[][][] tensor_estres;	

	// Variables relacionadas con la malla de nodos de fractura y los bordes
    boolean[][] 	vecinos_hook;
    int[]			borde_de_particula;
    List<Integer>	particulas_en_borde;
	double[][]  	posicion_nodos;
	double[][]  	vector_normal;
    double[][]  	vector_tangente;
    double[] 		punto_origen_fractura;
    double[]		punto_origen_fractura_recortado;
   	double[] 	 punto_borde_fractura;
	List<Integer> 	particulas_fractura;
	double[][]		traccion;
	double[]		superficie;
    List<double[]> 	fractura_origen;
    List<double[]> 	fractura_borde;
    
    
	//Otras variables
	double 		h_ini;
    double 		ancho_de_celda;
    double 		kernel_constante;
    double    	volumen;
    double[]  	separacion;
    double[]	angulo_estres;
    double[][]	estres_prin;
	double 		energia_elastica_total;
	double 		energia_cinetica_total;
    double[]	energia_auxiliar3;
  	double		trabajo_total;
	double[]	trabajo;
	double		energia_superficial_total;
	double[]	energia_superficial;
	double[]	tamano_inicial; 
	double[]	nuevos_tamanos;
	
	//Condiciones iniciales para el resto de las variables
	public  Condiciones_Iniciales()
	{
		
		//INICIACION DE LAS VARIABLES
		
		
        // Variables del metodo
        int N 					= Param.N;
        int tamano_inicial_vectores = 0;
        int t_actual 			= 0;
        double tiempo_anterior 	= Param.delta_t;
        double tiempo_inicial 	= 0; //primera iteracion del programa, se le especifica al metodo de tiempo que esta es la condicion inicial, es cero por definicion
        
        
        //Param de la transformacion Afin
        double[] vector_nulo 	= new double[dimension];
		double[] vector_unos 	= Util.Vector_Componentes_Iguales(dimension, 1.0);
        

        // Variables que transporta cada particula en equilibrio
        energia_elastica  		= new double[N];
        energia_cinetica  		= new double[N];
        trabajo					= new double[N];
		energia_superficial 	= new double[N];
        separacion    	 		= new double[dimension];
		velocidad     	 		= new double[N][dimension];
		densidad				= new double[N];
        
         // Variables que transporta cada particula relacionadas con las condiciobes de borde
        posicion      			= new double[N][dimension];
        posicion_sin_deformar 	= new double[tamano_inicial_vectores][dimension];
        posicion_nodos			= new double[N][2];
        elongacion    			= new double[N][dimension];
        tensor_deformacion 		= new double[N][dimension][dimension];
        tensor_estres 			= new double[N][dimension][dimension]; 
        
         //Variables relacionadas con el centro de masa  
        angulo_centro_de_masas 				= Param.angulo_centro_de_masas;
        posicion_particulas_centro_de_masa  = new double[tamano_inicial_vectores][dimension];
        posicion_centro_de_masa  			= Param.ubicacion;
	    velocidad_centro_de_masa 			= Velocidades_Iniciales(Param.velocidad_polares);
	    velocidad_particulas_centro_de_masa = new double[N][dimension];
        
         // Variables que transporta cada particula relacionadas con SPH
        h        						= new double[tamano_inicial_vectores];   
        vecinos  						= new int[N][]; 
        vecinos_hook  					= new boolean[N][]; 
        kernel   						= new double[N][];
        kernel_sin_consistencia   	 	= new double[N][];
        grad_kernel_sin_consistencia 	= new double[N][tamano_inicial_vectores][tamano_inicial_vectores]; 
        grad_kernel 					= new double[N][tamano_inicial_vectores][tamano_inicial_vectores]; 


		//Variables de las Tasa de cambio
        tasa_energia_elastica   				= new double[N];
        aceleracion   	 		= new double[N][dimension];

		// Otras variables, externas
		estres_prin  			= new double[N][dimension + 1];
		energia_auxiliar3		= new double[N];
		fractura_origen 		= new ArrayList<double[]>();
		fractura_borde 			= new ArrayList<double[]>();

		// Otras variables, internas
		double	 	 modulo_velocidad	= 0;
		double[] 	 tamano_fractura 	= Param.tamano_fractura;
		double[]  	 fuerza_externa 	= Otras_Var.Fuerza_Externa(dimension);
		double[]	 tamano 			= Param.tamano;
		double[] 	 h_inicial 			= new double[N];
        double[][] 	 velocidad_temporal = new double[N][dimension];
        
        
        //DEFINICION DE LAS VARIABLES



		Material Solido = new Material();
			volumen 							= Solido.volumen;
			tamano_inicial						= Solido.tamano_inicial;
			masas      							= Solido.masas;	
			energia_elastica 		    		= Solido.energia_elastica;						
			separacion							= Solido.separacion;
			posicion_particulas_centro_de_masa 	= Solido.posicion_particulas_centro_de_masa;				
			posicion_sin_deformar 				= Util.Ciclo_Vector_Trans_Afin(Solido.posicion_particulas_centro_de_masa, angulo_centro_de_masas, posicion_centro_de_masa, vector_nulo, vector_unos);
			velocidad_particulas_centro_de_masa = Solido.velocidad_particulas_centro_de_masa;

		//System.out.println(posicion_particulas_centro_de_masa[0][0] + " " + posicion_particulas_centro_de_masa[0][1]);
		//velocidad_particulas_centro_de_masa[12][0] = 0.1;
				
		numero_vecinos = (int)Math.pow(N, 1 - 2/Param.beta_escala_h);


		
	
		//Util.Imprimir_Vector(vecinos_hook[Param.particula_prueba]);

		Tipo_Condicion_Inicial Cond_Ini_Dinamica = new Tipo_Condicion_Inicial(masas, separacion, posicion_centro_de_masa, posicion_particulas_centro_de_masa, posicion_sin_deformar);
			posicion 	  				 = Cond_Ini_Dinamica.posicion;
			elongacion 	  				 = Cond_Ini_Dinamica.elongacion;
			kernel		  				 = Cond_Ini_Dinamica.kernel;
			grad_kernel   				 = Cond_Ini_Dinamica.grad_kernel;
			kernel_sin_consistencia 	 = Cond_Ini_Dinamica.kernel_sin_consistencia;
			grad_kernel_sin_consistencia = Cond_Ini_Dinamica.grad_kernel_sin_consistencia;
			tensor_estres 				 = Cond_Ini_Dinamica.tensor_estres;
			tensor_deformacion 			 = Cond_Ini_Dinamica.tensor_deformacion;
			densidad					 = Cond_Ini_Dinamica.densidad;
			vecinos						 = Cond_Ini_Dinamica.vecinos;
			vecinos_potenciales			 = Cond_Ini_Dinamica.vecinos_potenciales; 	
			celda_de_particula			 = Cond_Ini_Dinamica.celda_de_particula;
			ancho_de_celda				 = Cond_Ini_Dinamica.ancho_de_celda;
			h							 = Cond_Ini_Dinamica.h;


					
		long t_1 		= System.currentTimeMillis();

		/*
		 * mañana hacer la clase que haga lo siguiente
		 *  1) crear una linea de fractura
		 *  2) esta linea crea un rectangulo de ancho 2*separacion
		 * 
		 *  3) de las particulas que esten contenidas en este rectangulo, cuando se comparen
		 *  los vecinos, ver si la linea de fractura intercepta a un vecino, en caso de ser cierto
		 *  este vecino NO cuenta para el calculo del vector normal, esto se puede hacer en un ciclo
		 *  de j diferente al del ciclo j del vector normal de las particulas de borde, y el vector
		 *  normal resultante por particulas fractura se SUMA a las componentes del vector normal por
		 *  particulas borde, luego es que se normaliza, luego de eso es que se calcula el vector 
		 *  tangente y se añade esa particula a la lista de particulas_borde (teniendo cuidado que 
		 *  no se añada dos veces por borde y otra por fractura
		 * 
		 *  4) de este mismo ciclo i, se puede hacer la lista de vecinos_hook
		 * 
		 * 	5) cuando se crea una fractura, dos particulas se añaden a la lista del rectangulo/fractura 
		 *  y se calcula el vector normal/tangente y se añaden a la lista de particulas borde
		 * 
		 *  6) con todo esto, pensar como se haria el caso de particulas_borde en lo que se incluya una
		 *  cadena de particulas vecinas del borde incluyendo las fracturas (recordar que en el borde
		 *  de la fractura puede haber problemas
		 * 
		 *  7) para el calculo de la superficie, pensar bien como hacerlo con lo que tengo:
		 * 		- particulas del borde y sus vecinas
		 * 		- vector normal de cada una
		 *  con eso saber bien cuando identificar si la particula es de esquina o no, y con eso saber
		 *  si se se calcula la superficie o dos veces la superficie, tambien saber que con el area de 
		 *  un trapecio, conociendo la separacion entre los vecinos de una particula del borde y el
		 *  "area"/volumen de sph y la "profundidad"/separacion de la particula del borde de la particula
		 *  mas cercana del bulk, despejar la longitud faltante para el caso de una particula de la 
		 *  pared, para una de la esquina parece un poco mas complejo por que no es un trapecio
			*/
		
		
		double[] centro_factura_recta = Util.Transformacion_Afin_2D(Param.coordenadas_fractura, angulo_centro_de_masas, posicion_centro_de_masa, vector_nulo, vector_unos);
		
		Particulas_Fractura Par_Fac = new Particulas_Fractura(centro_factura_recta, tamano_fractura, posicion, separacion[0]); 

		punto_origen_fractura  			= Par_Fac.punto_origen_fractura;
		punto_borde_fractura			= Par_Fac.punto_borde_fractura;
		particulas_fractura 			= Par_Fac.particulas_fractura;
		posicion_nodos					= Par_Fac.rectangulo_particulas;
		//punto_origen_fractura_recortado = Par_Fac.punto_origen_fractura_recortado;
		//Util.Imprimir_Matriz(posicion_nodos);
	
		Particulas_En_Borde P_Borde_alt = new Particulas_En_Borde(separacion[0], h, vecinos, kernel_sin_consistencia, grad_kernel_sin_consistencia, particulas_fractura, punto_origen_fractura, punto_borde_fractura, posicion, masas, densidad);
		particulas_en_borde  	= P_Borde_alt.particulas_en_borde;
		vector_normal 			= P_Borde_alt.vector_normal;
		vector_tangente 		= P_Borde_alt.vector_tangente;
		vecinos_hook			= P_Borde_alt.vecinos_hook;
		//superficie				= P_Borde_alt.superficie;

		fractura_borde.add(punto_origen_fractura);
		fractura_borde.add(punto_borde_fractura);
		
		fractura_origen.add(punto_borde_fractura);
		fractura_origen.add(punto_origen_fractura);
		
		borde_de_particula		= Bordes.borde_de_particula(particulas_fractura, particulas_en_borde, separacion, posicion_particulas_centro_de_masa);
		traccion				= Condicion_Borde.Vector_Traccion(borde_de_particula, posicion_particulas_centro_de_masa);

		superficie = new Superficie_Particulas(particulas_en_borde, vecinos, posicion, vecinos_hook).superficie;
				
		for(int i =0; i<N; i++)
		{
			velocidad_temporal[i] 	= Util.Suma_Vectores(velocidad_particulas_centro_de_masa[i], velocidad_centro_de_masa);
			estres_prin[i] 			= Util.Autovalores_Reales_Matriz_2D(tensor_estres[i]);
			estres_prin[i][2]		= Otras_Var.Estres_Von_Mises(estres_prin[i][0], estres_prin[i][1]);	
		}
       		
       	Aceleracion Ace = new Aceleracion
		(	
			particulas_en_borde,
			borde_de_particula, 
			t_actual,
			vecinos, 
			masas, 
			densidad,
			superficie,
			grad_kernel,  
			tensor_estres,
			estres_prin, 
			vecinos_hook,
			kernel,
			h,
			traccion,
			vector_normal,
			posicion,
			separacion[0],
			velocidad_temporal,
			velocidad_temporal
		);	
			aceleracion = Ace.aceleracion;
		
		
		//Util.Imprimir_Matriz(estres_prin);
		//Util.Imprimir_Vector(borde_de_particula);
		
		tiempo  = new Tiempo(tiempo_anterior, h,  velocidad_temporal, aceleracion, vecinos).tiempo;
		//tiempo = 1;
		
		
		for(int i = 0; i< N; i++)	//calculo la fuerza corregida
		{
            velocidad[i]   		 	= Actual.Tensor_R1(velocidad_temporal[i], aceleracion[i], tiempo);
            modulo_velocidad 	 	= Util.Modulo_Vector(velocidad[i]);
			energia_cinetica[i]  	= 0.5*Param.densidad*Util.Producto_Punto(velocidad[i],velocidad[i]);
			energia_elastica[i]  	= 0.5*Util.Doble_Producto_Punto_Matrices(tensor_deformacion[i],tensor_estres[i] );
			nuevos_tamanos		 	= Otras_Var.Tamanos_Deformacion(tamano_inicial, tensor_deformacion[i]);
			trabajo[i]				= Otras_Var.Trabajo(masas, densidad[i], borde_de_particula[i], nuevos_tamanos, traccion[i],  elongacion[i], vector_normal[i]);
			energia_superficial[i] 	= Otras_Var.Energia_Superficial(masas, densidad[i], borde_de_particula[i], vector_normal[i], nuevos_tamanos);
		}
		
        tasa_energia_elastica   = Util.Tensor_Cero_Rango_1(N);
        tiempo  				= new Tiempo(tiempo, h,  velocidad, aceleracion, vecinos, 0).tiempo;
		//tiempo = 1;
		
		
		//velocidad[1][1] = 0.025;
		energia_elastica_total 		= Util.Suma_Componente_Vector(energia_elastica);
		energia_cinetica_total 		= Util.Suma_Componente_Vector(energia_cinetica);
		trabajo_total				= Util.Suma_Componente_Vector(trabajo);
		energia_superficial_total 	= Util.Suma_Componente_Vector(energia_superficial); 
		
	}
	

    public double[] Velocidades_Iniciales(double[] velocidad_polares)
    {
        int dimension  = Param.dimension;
        double PI_SOBRE_180 = 0.017453292;
        double angulo = PI_SOBRE_180*velocidad_polares[1];
        double modulo = Math.abs(velocidad_polares[0]);
        double[] velocidad = new double[dimension];
        velocidad[0] = modulo*Math.sin(angulo);
        velocidad[1] = modulo*Math.cos(angulo);
        return velocidad;
    }
    



    // revisa si la coordenada impuesta se ha salido de los bordes o esta en el borde
    public boolean[] Fuera_del_Borde(double[] posicion)
    {
        boolean[] fuera_del_borde = new boolean[Param.paredes + 1];
        fuera_del_borde[0] = posicion[0] <= Param.constante_izquierda;
        fuera_del_borde[1] = posicion[0] >= Param.constante_derecha;
        fuera_del_borde[2] = posicion[1] <= Param.constante_superior;
        fuera_del_borde[3] = posicion[1] >= Param.constante_inferior;
        fuera_del_borde[4] = fuera_del_borde[0] || fuera_del_borde[1] || fuera_del_borde[2] || fuera_del_borde[3];
        return fuera_del_borde;
    }

    // Regresa las particulas al dominio si por azar se generaron afuera de este
    public double[][] Posicion_Dentro_Borde(int N, double[][] posicion)
    {
        boolean[] condicion_borde = new boolean[Fuera_del_Borde(posicion[0]).length];
        double[][] posicion_corregida = new double[posicion.length][posicion[0].length];

        for(int i =0; i<N; i++)
        {
            condicion_borde = Fuera_del_Borde(posicion[i]);

            posicion_corregida[i][1] = posicion[i][1];
            posicion_corregida[i][0] = posicion[i][0];

            //evita que se salga de las paredes
            if (condicion_borde[4])
            {
                if ( condicion_borde[0] ) posicion_corregida[i][0] = Param.constante_izquierda;
                if ( condicion_borde[1] ) posicion_corregida[i][0] = Param.constante_derecha;
                if ( condicion_borde[2] ) posicion_corregida[i][1] = Param.constante_superior;
                if ( condicion_borde[3] ) posicion_corregida[i][1] = Param.constante_inferior;
            }
        }
        return posicion_corregida;
    }
    
	
}

class Tipo_Condicion_Inicial
{

    Utilidades 		Util 		= new Utilidades();
    Otras_Variables Otras_Var 	= new Otras_Variables();
    Parametros 		Param 		= new Parametros();
	Kernel_y_Derivadas Adap_y_Consis = new Kernel_y_Derivadas();
	
	double 		 ancho_de_celda;
	int[][] 	 vecinos;
	int[][] 	 vecinos_potenciales; 	
	int[][] 	 celda_de_particula;
	double[] 	 h;
	double[] 	 densidad;
	double[][] 	 posicion;
	double[][]	 elongacion;
	double[][]   kernel;
	double[][]   kernel_sin_consistencia;
	double[][][] grad_kernel;
	double[][][] grad_kernel_sin_consistencia;
	double[][][] tensor_deformacion;
	double[][][] tensor_estres;


	public Tipo_Condicion_Inicial(){}
	
	public Tipo_Condicion_Inicial
	(
		double 		masas,
		double[] 	separacion, 
		double[] 	posicion_centro_de_masa,
		double[][] 	posicion_particulas_centro_de_masa, 
		double[][]  posicion_sin_deformar
	)
	{

		int 	 	N 						= posicion_sin_deformar.length;
		int 	 	dimension 				= posicion_sin_deformar[0].length;
		int 	 	tipo_cond_inicial 		= Param.tipo_cond_inicial;
		int 	 	tamano_inicial_vectores = 0;
		int[][]  	vecinos_ini			 	= new int[N][1];
		double[]	tamano					= Param.tamano;
		double 		h_ini 					= H_Inicial( N, separacion, tamano);
        double[]    h_inicial				= Util.Vector_Componentes_Iguales(N, h_ini);
                
        h							 = new double[N];
		densidad					 = new double[N];
		posicion 					 = new double[N][dimension];
		elongacion  				 = new double[N][dimension];
		tensor_deformacion 			 = new double[N][dimension][dimension];
		tensor_estres 	   			 = new double[N][dimension][dimension];   
		kernel   					 = new double[N][tamano_inicial_vectores];
        grad_kernel 				 = new double[N][tamano_inicial_vectores][tamano_inicial_vectores]; 		
        kernel_sin_consistencia      = new double[N][tamano_inicial_vectores];    
        grad_kernel_sin_consistencia = new double[N][tamano_inicial_vectores][tamano_inicial_vectores];  
		
		

		
		//System.out.println(posicion[0][0] + " " + posicion_sin_deformar[0][0]);
				
		for(int i =0; i < N; i++)
		{
		//	elongacion[i] = Util.Resta_Vectores(posicion[i], posicion_sin_deformar[i]);
			
			tensor_estres[i] 	  = Param.estres_inicial;
			tensor_deformacion[i] = Otras_Var.Tensor_Deformacion_A_Partir_Estres(i, tensor_estres[i]);
			elongacion[i] 		  = Otras_Var.Elongaciones_Deformacion_Constante(separacion[0], new double[dimension], posicion_particulas_centro_de_masa[i],  tensor_deformacion[i]);
			posicion[i] 		  = Util.Suma_Vectores(elongacion[i], posicion_sin_deformar[i]);
			densidad[i]			  = Param.densidad;	//el vector densidad usado para el calculo de sph es una densidad numerica por unidad de area, la densidad dada en parametros es la densidad real del material por unidad de volumen

		}
		

		//System.out.println(tensor_deformacion[0][0][0] + " " +  tensor_deformacion[Param.particula_prueba][0][0]);
		
		Busqueda_Vecinos Vecinos_Prueba = new Busqueda_Vecinos(h_inicial, posicion);
		vecinos_ini        = Vecinos_Prueba.vecinos_reales;
		
		//System.out.println(h_ini);
		
		h = Adap_y_Consis.H(posicion, masas, densidad, h_inicial, vecinos_ini);
		
		//System.out.println(h[0]);
		
		Busqueda_Vecinos Vecinos = new Busqueda_Vecinos(h, posicion);
			vecinos        		= Vecinos.vecinos_reales;
			ancho_de_celda 		= Vecinos.ancho_de_celda;
			vecinos_potenciales = Vecinos.vecinos_potenciales; 	
			celda_de_particula 	= Vecinos.celda_de_particula; 	
			

		for(int i =0; i<N; i++)
		{
			Kernel_y_Derivadas Kernel 			 = new Kernel_y_Derivadas(i, h, vecinos[i], posicion);
				kernel_sin_consistencia[i]  	 = Kernel.kernel;
				grad_kernel_sin_consistencia[i]  = Kernel.gradiente;
		}
		
		
		for(int i = 0; i< N; i++) // kernel y sus derivadas con consistencia orden 0
		{
			kernel[i] 	   = Adap_y_Consis.Cons_Kernel(i, vecinos[i], masas, densidad, kernel_sin_consistencia[i], posicion);
			grad_kernel[i] = Adap_y_Consis.Cons_Grad(i, vecinos[i], masas, densidad, posicion, kernel_sin_consistencia[i], kernel[i], grad_kernel_sin_consistencia[i]);
		}
		
	}
	
	 public double Kernel_Constante(double h, double separacion, int dimension)
    {
        double q = separacion/(h);
        return  new Kernel_y_Derivadas().Kernel(-1, -1, q, h, dimension);
    }
    
    public double H_Inicial(int N_bloque, double[] separacion, double[] tamano)
    {
		boolean h_separacion 			= Param.h_separacion;
        double  beta_escala_h 			= Param.beta_escala_h;
        double  constante_num_vecinos	= 4;
        
		if(h_separacion)
			return Param.parametro_D*separacion[0];
		else
		{
			//System.out.println(tamano[0] + " " + tamano[1] + " " + N_bloque  + " " + Math.sqrt(tamano[0]*tamano[1]/Math.PI)*Math.pow(N_bloque,-1/beta_escala_h));
			if(N_bloque > 2)
				return Math.sqrt(tamano[0]*tamano[1]*constante_num_vecinos)/Math.PI*Param.h0;
				//return 1.5*separacion[0];
			else
				return 2*Math.sqrt(tamano[0]*tamano[1]/Math.PI)*Param.h0;
		}
    }
}



