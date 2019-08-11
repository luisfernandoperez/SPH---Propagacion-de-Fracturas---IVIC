import java.io.*;
import java.util.*;
import java.util.List;
import java.awt.*;
import java.lang.Runtime;
import javax.swing.JOptionPane;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.Timer;

class Main_SPH extends JPanel implements ActionListener
{
	//objetos de los metodos
	static Utilidades 	Util  = new Utilidades();
    static Parametros	Param = new Parametros();

	//variables de los tiempos de animacion y duracion del programa
	static int 	tmax 			= Param.tmax;
	static long t_inicial 		= System.currentTimeMillis();
	static int 	t_actual 		= 0;
	static double tiempo_real;
	double  	tiempo_acustico = Param.tiempo_acustico;
	int    		ms_entre_frames = Param.ms_entre_frames;
	Timer  		tm 				= new Timer(ms_entre_frames, this);
	int 		t_anterior 		= -1;
	double 		delta_tiempo_inicial;
	double 		delta_tiempo;
	double 		delta_tiempo_anterior;
	boolean 	tiempo_max_acus = Param.tiempo_max_acus;


	//Variables propias de SPH
	double[]	 h;
	double	 	 h_inicial;
	int   [][]	 vecinos;
	int[][] 	 vecinos_potenciales; 	
	int[][] 	 celda_de_particula;
    double[][]	 kernel;
    double[][][] grad_kernel;


	//Variables que se actualizan en cada paso de tiempo
	double[] 	 energia_elastica ;
	double[]	 energia_cinetica;
	double[]	 densidad;
	double[][]	 posicion;
	double[][]	 velocidad;
	double[][][] tensor_estres;

	//Tasas de cambio de cada una de las variables de arriba
	double[] 	 tasa_energia_elastica ;
	double[][]	 aceleracion;

	//Variables predichas o intermedias del leapfrog

	double[]	 tasa_energia_elastica_predicha;
	double[][]   velocidad_intermedia;
			
	//Variables anteriores, variables temporales
	double[] 	 energia_elastica_anterior;
	double[][]	 posicion_anterior;
	double[][]	 velocidad_anterior;
	
	//Variables predichas
	double[]	 energia_elastica_predicha;

    //Tasas de cambio anteriores, variables temporales
    double[]     tasa_energia_elastica_anterior;
    double[][]   aceleracion_anterior;
    
    //Variables relacionadas con el centro de masa
    double[] 	posicion_centro_de_masa;
    double[]	posicion_centro_de_masa_inicial;
	double[][] 	posicion_particulas_centro_de_masa;
	double[] 	velocidad_centro_de_masa;
	double[][]	velocidad_particulas_centro_de_masa;
    
    //Variables relacionadas con las energias
	double 		energia_cinetica_total;
	double 		energia_elastica_total;
	double[] 	energia_auxiliar3;
	double		trabajo_total;
	double[]	trabajo;
	double		energia_superficial_total;
	double[]	energia_superficial;
		
	//Otras Variables fisicas 
	double 	 	 masas;
	double 	 	 volumen;
    double[][]   posicion_inicial;
	double[][] 	 estres_prin;
    double[][] 	 elongacion;
    double[][][] tensor_deformacion;
	
	
	 //Variables relacionadas con los bordes y las fracturas
    List<Integer>	particulas_en_borde;
    List<Integer>	particulas_fractura;
	boolean[][]  	vecinos_hook;
	int[]		 	borde_de_particula;
	double[][]   	vector_normal;
	double[][]   	vector_tangente;
    double[] 		punto_origen_fractura;
    double[] 		punto_borde_fractura;
    double[][]		traccion;
    double[]		separacion;
	double[][]		posicion_nodos;
	double[]		superficie;
	double			longitud_fractura;
	double			longitud_fractura_anterior;
	double[] 		tamano_inicial;
	double			criterio_de_griffith;
	List<double[]> 	fractura_origen;
	List<double[]> 	fractura_borde;
	boolean			llego_origen_a_frontera;
	boolean			llego_borde_a_frontera;
		
		
	//Otras variables
	int 		N = Param.N;
	int		    numero_vecinos;
	int 		particula_prueba = Param.particula_prueba;
	double 	 	r;
	double		ancho_de_celda;
	boolean solo_una_vez = true;
					
    public static void main(String args[])
    {
		System.out.println(" ");	//imprime una linea blanca solo por estetica
		JFrame ventana = new JFrame(Param.titulo);
        ventana.setResizable(false);
        ventana.getContentPane().add(new Main_SPH());
        ventana.setSize(Param.x_lienzo ,Param.y_lienzo);
        ventana.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        ventana.setLocationRelativeTo(null);
		ventana.addWindowListener
		(
			new java.awt.event.WindowAdapter() 
			{
				@Override
				public void windowClosing(java.awt.event.WindowEvent windowEvent) 
				{				
					Util.Limpiar_Archivos(new File(Param.directorio_base), "class");
					if(tmax< 0 || t_actual < tmax)
						Finalizacion(t_inicial, t_actual);
					System.exit(0);
				}
			}
		);
       
        ventana.setVisible(true);
    }
      
    @Override
    public void paintComponent(Graphics g)
    {

		Inicializacion();

		if(t_actual != t_anterior )
		{ 
			Colocar_Condiciones_Iniciales(t_actual);
			Ciclo_Tiempo(t_actual);
			
			new Imprimir
			(
				t_actual,
				tmax,
				tiempo_real,
				posicion,
				elongacion,	  	  
				posicion_centro_de_masa,  	  
				posicion_particulas_centro_de_masa,
				velocidad,
				velocidad_centro_de_masa,	  
				velocidad_particulas_centro_de_masa,	    
				aceleracion,  
				energia_elastica ,
				tasa_energia_elastica ,		   
				tensor_deformacion,	  
				tensor_estres,  	  
				estres_prin,	  		
				h,				
				energia_cinetica_total,
				energia_elastica_total,
				energia_auxiliar3,
				vecinos_hook,
				densidad,
				punto_origen_fractura,
				separacion[0],
				trabajo_total,
				trabajo,
				energia_superficial_total,
				energia_superficial,
				energia_cinetica,
				criterio_de_griffith,
				longitud_fractura,
				longitud_fractura_anterior
		);

			if(Param.graficacion)
			{
				super.paintComponent(g);
				Graphics2D g_2D = (Graphics2D) g; 
				
				Graficar graf = new Graficar
				(
					g_2D,
					t_actual,
					tmax,
					tiempo_real,
					posicion,
					elongacion,	  	  
					posicion_centro_de_masa,  	  
					posicion_particulas_centro_de_masa,
					velocidad,
					velocidad_centro_de_masa,	  
					velocidad_particulas_centro_de_masa,	  
					aceleracion,  
					energia_elastica ,
					tasa_energia_elastica ,		  
					tensor_deformacion,	  
					tensor_estres,  
					estres_prin,	  		
					h,
					ancho_de_celda,
					vecinos,
					particulas_en_borde,
					r,
					vecinos_hook,
					posicion_nodos,
					densidad,
					borde_de_particula,
					fractura_origen,
					fractura_borde
				);
				
					graf.Lienzo();
					graf.Auto_Ajuste();
					graf.Mallado();
					graf.Barra_Lateral();
					graf.Graficar_Particulas();
					graf.Graficar_Particulas_En_Especifico();
					graf.Graficar_Nodos();
					graf.Graficar_Numeros_Particulas();
					graf.Graficar_Fractura();
			}
			
			if(t_actual == tmax)	//finalizacion del programa
			{	
				Finalizacion(t_inicial, t_actual);
			}
		}
		t_anterior = t_actual;
	}

    
    @Override		//animacion
    public void actionPerformed(ActionEvent e)
    {
		// termina el ciclo del tiempo, si tmax es negativo, la animacion no termina
        if (t_actual >= tmax  && tmax > -1)	
        {
            tm.stop();
            return;
        }
        else
        {
			t_actual++;
            repaint();
        }
    }

	//Inicia la animacion
    public void Inicializacion()
    {
        tm.start();    
	}
	
	public void Colocar_Condiciones_Iniciales(int t_actual)
	{ 
		if(t_actual == 0)	//Condiciones iniciales
		{

			//Condiciones iniciales de las  coordenadas, sus velocidades y las variables de estado
			Condiciones_Iniciales condiciones = new Condiciones_Iniciales();

				//Condiciones iniciales de variables basicas
				masas  					= condiciones.masas;
				volumen 				= condiciones.volumen;
				densidad				= condiciones.densidad;
				delta_tiempo 			= condiciones.tiempo;
				delta_tiempo_inicial	= delta_tiempo;
				tamano_inicial			= condiciones.tamano_inicial;
								
				//Condiciones iniciales de variables de la integracion mecanica
				posicion 			= condiciones.posicion;
				velocidad 			= condiciones.velocidad;
				aceleracion 		= condiciones.aceleracion;
				tensor_estres 		= condiciones.tensor_estres;
				tensor_deformacion 	= condiciones.tensor_deformacion;
				elongacion 			= condiciones.elongacion;
										
				//Condiciones iniciales relacionadas intrisecamente con SPH
				h	  				= condiciones.h;
				h_inicial 			= condiciones.h_ini;
				vecinos 			= condiciones.vecinos;
				vecinos_potenciales	= condiciones.vecinos_potenciales;
				celda_de_particula	= condiciones.celda_de_particula;
				kernel 				= condiciones.kernel;
				grad_kernel 		= condiciones.grad_kernel;
				

				//Variables relacionadas con el centro de masa
				posicion_particulas_centro_de_masa  = condiciones.posicion_particulas_centro_de_masa;   
				posicion_centro_de_masa 			= condiciones.posicion_centro_de_masa;
				posicion_centro_de_masa_inicial		= condiciones.posicion_centro_de_masa;
				velocidad_centro_de_masa 			= condiciones.velocidad_centro_de_masa;
				velocidad_particulas_centro_de_masa = condiciones.velocidad_particulas_centro_de_masa;
				posicion_inicial 					= condiciones.posicion_sin_deformar;
				 
				//Variables relacionadas con los bordes y las fracturas 
				particulas_en_borde		= condiciones.particulas_en_borde;
				particulas_fractura     = condiciones.particulas_fractura;
				borde_de_particula 		= condiciones.borde_de_particula;			
				vector_normal 			= condiciones.vector_normal;			
				vector_tangente 		= condiciones.vector_tangente;			
				punto_origen_fractura	= condiciones.punto_origen_fractura;
				punto_borde_fractura	= condiciones.punto_borde_fractura;
				traccion				= condiciones.traccion;
				posicion_nodos			= condiciones.posicion_nodos;
				superficie				= condiciones.superficie;
				
				//condicion iniciales relacionadas con la energia
				energia_elastica  		= condiciones.energia_elastica;
				energia_cinetica 				= condiciones.energia_cinetica;
				trabajo						= condiciones.trabajo;
				energia_elastica_total 		= condiciones.energia_elastica_total;
				energia_cinetica_total 		= condiciones.energia_cinetica_total;
				trabajo_total				= condiciones.trabajo_total;
				energia_superficial_total 	= condiciones.energia_superficial_total;
				energia_superficial			= condiciones.energia_superficial;
				energia_auxiliar3 			= condiciones.energia_auxiliar3;
				tasa_energia_elastica  		= condiciones.tasa_energia_elastica; //por ahora no implementado
				
				//Condiciones iniciales de variables auxiliares
				separacion 		 		= condiciones.separacion;
				ancho_de_celda 	 		= condiciones.ancho_de_celda*Param.escala_lienzo;
				vecinos_hook  			= condiciones.vecinos_hook;
				r 						= condiciones.separacion[0];
		        estres_prin				= condiciones.estres_prin;
		        numero_vecinos			= condiciones.numero_vecinos;
		        fractura_origen			= condiciones.fractura_origen;
		        fractura_borde			= condiciones.fractura_borde;
		        
				longitud_fractura		= Param.tamano_fractura[0];
				longitud_fractura_anterior = longitud_fractura;
				llego_origen_a_frontera = false;
				llego_borde_a_frontera 	= false;
				criterio_de_griffith 	= new Criterio_Griffith(longitud_fractura).criterio_de_griffith;
				//Util.Imprimir_Lista(particulas_fractura);

				System.out.println(Param.vel_sonido*Math.sqrt(Param.modulo_compresibi_real/Param.densidad_real));
				
				
		}
	}
	

	public void Ciclo_Tiempo(int t_actual)
	{
		if(t_actual != 0  && (t_actual <= tmax || tmax < 0))	//leapfrog
		{
			
			tiempo_real += delta_tiempo;
			
			if(tiempo_max_acus && tiempo_real > tiempo_acustico)
			{
				Finalizacion(t_inicial, t_actual);
				tm.stop();
			}
			
            //Variables temporales que guardan los valores de las variables antes de ser modificados
            longitud_fractura_anterior		= longitud_fractura;
			delta_tiempo_anterior		   	= delta_tiempo;
			energia_elastica_anterior 		= energia_elastica ;
			posicion_anterior 	   			= posicion;
			velocidad_anterior 	   			= velocidad;
            tasa_energia_elastica_anterior 	= tasa_energia_elastica ;
            aceleracion_anterior   	  		= aceleracion;

			//System.out.println(Param.criterio_de_griffith);
			
			//Integrador numerico: primera etapa
			Verlet Verlet_Inicial = new Verlet
			(
				delta_tiempo_anterior,
				particulas_en_borde, 
				borde_de_particula,
				posicion_anterior,
				posicion_inicial, 
				posicion_centro_de_masa_inicial,
				velocidad_anterior, 
				aceleracion_anterior
			);
													 
				velocidad_intermedia 				= Verlet_Inicial.velocidad_intermedia;
				posicion   							= Verlet_Inicial.posicion;
				elongacion 							= Verlet_Inicial.elongacion;
				posicion_centro_de_masa 			= Verlet_Inicial.posicion_centro_de_masa;
				posicion_particulas_centro_de_masa 	= Verlet_Inicial.posicion_particulas_centro_de_masa;

			
			//posicion[0][1] = posicion_anterior[0][1];
			//posicion[2][1] = posicion_anterior[2][1];

			//Integrador numerico
            Corrector Paso_Corrector = new Corrector
            (
				t_actual,
				delta_tiempo_anterior,
				tiempo_real, 
				masas, 
				energia_elastica_anterior, 
				tasa_energia_elastica_anterior,
				tasa_energia_elastica_predicha,	
				posicion, 
				elongacion, 
				velocidad_intermedia, 
				h, 
				vecinos,
				vecinos_potenciales, 
				celda_de_particula,
				particulas_en_borde,
				borde_de_particula, 
				aceleracion_anterior,
				separacion,
				vecinos_hook,
				traccion,
				densidad,
				particulas_fractura,
				punto_origen_fractura,
				vector_normal,
				tamano_inicial,
				posicion_particulas_centro_de_masa,
				superficie,
				punto_borde_fractura,
				longitud_fractura,
				velocidad_anterior,
				fractura_borde,
				fractura_origen,
				llego_origen_a_frontera,
				llego_borde_a_frontera
			);

				
				delta_tiempo 	    		= Paso_Corrector.delta_tiempo;
				velocidad 					= Paso_Corrector.velocidad;
				aceleracion 				= Paso_Corrector.aceleracion;
				tensor_estres 				= Paso_Corrector.tensor_estres;
				tensor_deformacion 			= Paso_Corrector.tensor_deformacion;
				
				
				estres_prin					= Paso_Corrector.estres_prin;
				
				energia_elastica	 		= Paso_Corrector.energia_elastica;
				energia_cinetica 			= Paso_Corrector.energia_cinetica;
				trabajo						= Paso_Corrector.trabajo;
				energia_superficial			= Paso_Corrector.energia_superficial;
				
				energia_elastica_total 		= Paso_Corrector.energia_elastica_total;
				energia_cinetica_total 		= Paso_Corrector.energia_cinetica_total;
				trabajo_total				= Paso_Corrector.trabajo_total;
				energia_superficial_total 	= Paso_Corrector.energia_superficial_total;
				
	
				vecinos_hook 				= Paso_Corrector.vecinos_hook_nuevo; 
				punto_origen_fractura		= Paso_Corrector.punto_origen_fractura_nuevo;
				punto_borde_fractura		= Paso_Corrector.punto_borde_fractura_nuevo;
				longitud_fractura			= Paso_Corrector.longitud_fractura_nuevo;
				criterio_de_griffith		= Paso_Corrector.criterio_de_griffith;
				fractura_borde				= Paso_Corrector.fractura_borde_nuevo;
				fractura_origen				= Paso_Corrector.fractura_origen_nuevo;
				particulas_en_borde			= Paso_Corrector.particulas_en_borde_nuevo;
				vector_normal				= Paso_Corrector.vector_normal_nuevo;
				borde_de_particula			= Paso_Corrector.borde_de_particula_nuevo;
				superficie					= Paso_Corrector.superficie_nuevo;
				llego_origen_a_frontera		= Paso_Corrector.llego_origen_a_frontera_nuevo;
				llego_borde_a_frontera		= Paso_Corrector.llego_borde_a_frontera_nuevo;
		
				//Util.Imprimir_Vector(punto_origen_fractura);
				
			velocidad_centro_de_masa = Util.Vector_Centro_de_Masas(velocidad);
			
				//Util.Imprimir_Lista(particulas_fractura);
						
		}
	}





	static public void Finalizacion(long t_inicial, int t_actual)
	{

		long t_final = System.currentTimeMillis();
		float segundos = (t_final - t_inicial)/1000f;
		int minutos = (int)(segundos/60f);
		if(segundos>= 60.) segundos -= minutos*60f;
		int horas = (int)(minutos/60f);
		if(minutos >= 60) minutos -= horas*60;
		int dias = horas/24;
		if(horas >= 24) horas -= dias*24;
		System.out.println(" ");
		System.out.println("El programa ha finalizado en la iteracion numero " + (t_actual) + " en");
		System.out.println(dias + " dias " + horas + " horas " + minutos + " minutos y " + segundos + " segundos" );
		System.out.println("y el tiempo de simulacion ha sido " + tiempo_real*Param.escala_tiempo_real + " segundos");
		System.out.println(" ");
		System.out.println(" ");
		Util.Limpiar_Archivos(new File(Param.directorio_base), "class");
		if(Param.cerrar_despues_finalizar)
			System.exit(0);
	}
	
	



    
}
