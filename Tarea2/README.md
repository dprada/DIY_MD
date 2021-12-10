# Tarea 2

## Objetivo

Programar la integración de la dinámica de Langevin de una partícula en tres dimensiones en un
potencial externo cualquiera -en su versión más simple y erronea (Euler)-:

\begin{align}
\vec{r}(t_{i+1}) &= \vec{r}(t_{i}) + \vec{v}(t_{i})\Delta t \\
\vec{v}(t_{i+1}) &= \vec{v}(t_{i}) + \vec{a}(t_{i})\Delta t
\end{align}

donde:

\begin{equation}
m\vec{a}(t_{i}) = \sum_{l}{\vec{F_l}(\vec{r}(t_{i}), \vec{v}(t_{i}))}
\end{equation}

## Proceso

### Paso 1 (Dinámica clásica)

Programar la integración de la dinámica clásica de un partícula de masa m en un potencial externo
cualquiera. El resultado debe ser la trayectoria, posiciones y velocidades, en función del tiempo
almacenada en memoria como dos numpy arrays.

El potencial, o las fuerzas externas, deben estar definidas en su propia función de Python.

El paso de integración, la masa, la posición inicial, las velocidades iniciales y el número de
pasos a integrar (o el tiempo total de simulación) deben ser parámetros libres del problema
(variables que poder modificar). 

Como propuesta de potencial externo inicial optaremos por:

\begin{equation}
V(x,y,z) = \frac{1}{2} k_{x} x^{2} + \frac{1}{2} k_{y} y^{2} + \frac{1}{2} k_{z} z^{2}
\end{equation}

### Paso 2 (Dinámica clásica amortiguada)

Al paso anterior le vamos a incluir una fuerza externa de amortiguamiento dependiente de un
parametro gamma -damping- de dimensión física 1/T (T: tiempo). La fuerza por tanto tiene la forma:

\begin{equation}
F_{damping} = - \gamma m\vec{v}
\end{equation}

El parámetro $\gamma$ será por tanto una variable libre de entrada de nuestra rutina. Para poder si
modificada y comprobar el resultado de la simulación.

### Paso 3 (Dinámica de Langevin)

Vamos a introducir en el resultado del paso dos el término de fluctuación de la dinámica de
Langevin. Para ello necesitaremos incluir el término $K_{B}T$ como variable en la rutina, donde $K_{B}$ es la constante de Boltzmann y $T$ es la temperatura (en Kelvin).


