import java.io.*;

public class Tiempo
{
    Utilidades Util = new Utilidades();
    Parametros Param = new Parametros();
	double tolerancia = 0.01;
	double tiempo;

	public Tiempo
	(
		double tiempo_anterior, 
		double[] h, 
		double[][] velocidad, 
		double[][] aceleracion, 
		int[][] vecinos, 
		int t_actual
	)
	{
	boolean tiempo_variable	= Param.tiempo_variable;
	int N  = Param.N;
	double tiempo_fuerza = 0;
	double tiempo_cfl = 0;
	double tiempo_acustico = 0;
	double tiempo_cf = 0;
	double t_fue_temporal = 0;
	double t_cfl_temporal = 0;
	double t_acu_temporal = 0;
	double t_cf_temporal = 0;
	double cons_tiempo = Param.cons_tiempo;
	double t_temporal_1 = 0;
	double t_temporal_2 = 0;
    double const1 = 0.2;
    double const2 = 0.1;
    double tiempo_minimo = 0.0; // esta constante existe para que este metodo calcule el valor minimo de dos tiempos diferentes pero que sea distinto de cero
	boolean condicion_cfl;
	boolean condicion_fuerza;
	boolean condicion_cf;
	int n = 1159;

//System.out.println( h[n] + " " + velocidad[n][0] + " " + velocidad[n][1] + " " + t_actual);
//Util.Imprimir_Matriz(velocidad);

		if(tiempo_variable)	//Este ciclo es para calcular el tiempo, por optimizacion se puede incluir en otros ciclos
		{
			tiempo_fuerza	= Tiempo_Fuerza(h[0], aceleracion[0]);
			tiempo_cfl		= Tiempo_CFL(t_actual, 0, h[0], velocidad, vecinos[0]);
			tiempo_acustico = Tiempo_Acustico(h[0], Param.vel_sonido);
			tiempo_cf		= Tiempo_CF(h[0], velocidad[0]);
			
			for(int i = 1; i< N; i++)
			{
				t_fue_temporal	= Tiempo_Fuerza(h[i], aceleracion[i]);
				t_cfl_temporal	= Tiempo_CFL(t_actual, i, h[i], velocidad, vecinos[i]);
				t_acu_temporal	= Tiempo_Acustico(h[i], Param.vel_sonido);
				t_cf_temporal	= Tiempo_CF(h[i], velocidad[i]);
				//if(i == n)System.out.println(t_fue_temporal + " " + t_cfl_temporal + " " + t_acu_temporal + " " + t_cf_temporal);
				tiempo_fuerza	= Util.Maximo_2_numeros(tiempo_minimo, Util.Minimo_2_numeros(tiempo_fuerza, t_fue_temporal));
				tiempo_cfl		= Util.Maximo_2_numeros(tiempo_minimo, Util.Minimo_2_numeros(tiempo_cfl, t_cfl_temporal));
				tiempo_acustico = Util.Maximo_2_numeros(tiempo_minimo, Util.Minimo_2_numeros(tiempo_acustico, t_acu_temporal));
				tiempo_cf		= Util.Maximo_2_numeros(tiempo_minimo, Util.Minimo_2_numeros(tiempo_cf, t_cf_temporal));
				
				//if(i == n)System.out.println( h[n] + " " + velocidad[n][0] + " " + velocidad[n][1] + " " + aceleracion[n][0] + " " + aceleracion[n][1]);
				condicion_cfl 	 = tiempo_cfl != tiempo_cfl;
				condicion_fuerza = tiempo_fuerza != tiempo_fuerza;
				condicion_cf 	 = tiempo_cf != tiempo_cf;
				
				if(condicion_cfl || condicion_fuerza || condicion_cf)
				{ 
					
					if(condicion_cfl)
						System.out.println("hay un problema en el metodo de tiempos para la iteracion " + t_actual + " con la particula " + i + " para el metodo de tiempo CFL");
					if(condicion_fuerza)
						System.out.println("hay un problema en el metodo de tiempos para la iteracion " + t_actual + " con la particula " + i + " para el metodo de tiempo de fuerzas");
					if(condicion_cf)
						System.out.println("hay un problema en el metodo de tiempos para la iteracion " + t_actual + " con la particula " + i + " para el metodo de tiempo CF");
					
					Util.Limpiar_Archivos(new File(Param.directorio_base), "class");
					System.exit(1);
				}
			}

			t_temporal_1 = Util.Minimo_2_numeros(const1*tiempo_fuerza, tiempo_acustico);
			t_temporal_2 = Util.Minimo_2_numeros(tiempo_cf, tiempo_cfl);
			tiempo = cons_tiempo*Util.Minimo_2_numeros(t_temporal_1, const2*t_temporal_2);
			//if(tiempo < Param.delta_t) tiempo = Param.delta_t;
			if(tiempo != tiempo) tiempo = tiempo_anterior;
			//if(tiempo < 0.00001) tiempo =0.00001;

		}
		if(!tiempo_variable){tiempo = tiempo_anterior;} 
		//System.out.println(tiempo + " " + (float)tiempo_fuerza + " " + (float)tiempo_cfl + " " + (float)tiempo_acustico + " " + (float)tiempo_cf );
    }

    public double Tiempo_Fuerza(double h, double[] aceleracion)
    {
		//if(h < h_min) h = h_min;
        double modulo =  Util.Modulo_Vector(aceleracion);
		if(modulo < tolerancia) modulo = tolerancia;
        return Math.sqrt(h/modulo);
    }

	 public double Tiempo_CFL(int t_actual, int i, double h, double[][] velocidad, int[] vecinos)
    {
		//if(h < h_min) h = h_min;
        double tiempo_adaptativo = 0;
		double[] v_ij = new double[velocidad[0].length];
		double modulo = 0;
		double temporal = 0;
		//if(t_actual == 1 && i == 311) Util.Imprimir_Vector(velocidad[311]);
        for(int j = 0; j<vecinos.length; j++)
		{	
			v_ij = Util.Resta_Vectores(velocidad[i], velocidad[vecinos[j]]);
			modulo = Util.Modulo_Vector(v_ij);
			if(modulo < tolerancia) modulo = tolerancia;
			//if(t_actual == 1 && i == 915) System.out.println(vecinos[j] + " " + modulo + " " + velocidad[vecinos[j]][0] + " " + velocidad[vecinos[j]][1]);
			temporal = h/modulo;
			if(j== 0) tiempo_adaptativo = temporal;
			else tiempo_adaptativo = Util.Minimo_2_numeros(tiempo_adaptativo, temporal);
		}
        return tiempo_adaptativo;
    }

	 public double Tiempo_CF(double h, double[] velocidad)
    {
		//if(h < h_min) h = h_min;
		double modulo = 0;
		modulo = Util.Modulo_Vector(velocidad);
		if(modulo < tolerancia) modulo = tolerancia;
        return h/modulo;
    }

	 public double Tiempo_Acustico(double h, double velocidad_sonido)
    {
		//if(h < h_min) h = h_min;
		if(velocidad_sonido < tolerancia) velocidad_sonido = tolerancia;
        return h/(velocidad_sonido*(1+ 0.0*Param.alfa_viscocidad_artificial));
    }

        //este solo se usa en las condiciones iniciales
    public Tiempo(double tiempo_anterior, double[] h, double[][] velocidad, double[][] aceleracion, int[][] vecinos)
    {
        boolean tiempo_variable    = Param.tiempo_variable;
        int N  = Param.N;
        double tiempo_fuerza;
        double tiempo_cfl;
        double tiempo_acustico;
        double tiempo_cf;
        double t_fue_temporal;
        double t_cfl_temporal;
        double t_acu_temporal;
        double t_cf_temporal;
        double cons_tiempo = Param.cons_tiempo;
        double t_temporal_1;
        double t_temporal_2;

        tiempo_fuerza    = Tiempo_Fuerza(h[0], aceleracion[0]);
        tiempo_cfl        = Tiempo_CFL(0, 0, h[0], velocidad, vecinos[0]);
        tiempo_acustico = Tiempo_Acustico(h[0], Param.vel_sonido);
        tiempo_cf        = Tiempo_CF(h[0], velocidad[0]);
        for(int i = 1; i< N; i++)
        {
            t_fue_temporal  = Tiempo_Fuerza(h[i], aceleracion[i]);
            t_cfl_temporal  = Tiempo_CFL(0, i, h[i], velocidad, vecinos[i]);
            t_acu_temporal  = Tiempo_Acustico(h[i], Param.vel_sonido);
            t_cf_temporal   = Tiempo_CF(h[i], velocidad[i]);
            tiempo_fuerza   = Util.Minimo_2_numeros(tiempo_fuerza, t_fue_temporal);
            tiempo_cfl      = Util.Minimo_2_numeros(tiempo_cfl, t_cfl_temporal);
            tiempo_acustico = Util.Minimo_2_numeros(tiempo_acustico, t_acu_temporal);
            tiempo_cf       = Util.Minimo_2_numeros(tiempo_cf, t_cf_temporal);
        }
        t_temporal_1 = Util.Minimo_2_numeros(0.3*tiempo_fuerza, tiempo_acustico);
        t_temporal_2 = Util.Minimo_2_numeros(tiempo_cf, tiempo_cfl);
        tiempo = cons_tiempo*Util.Minimo_2_numeros(t_temporal_1, t_temporal_2);
        //if(tiempo < Param.delta_t) tiempo = Param.delta_t;
        if(tiempo != tiempo) tiempo = Param.delta_t;
        //if(tiempo < 0.00001) tiempo =0.00001;
       // System.out.println( (float)tiempo_fuerza + " " + (float)tiempo_cfl + " " + (float)tiempo_acustico + " " + (float)tiempo_cf );

    }
}

