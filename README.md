# C-basic

Recopilación de códigos sencillos escritos en C.
Todos los códigos están comentados con el objetivo de servir como material de aprendizaje.

> Cada sección está plegada. Haz clic en un título para desplegarlo.

---

<details open>
<summary><b>🚀 Cómo compilar y ejecutar</b></summary>

<br>

Para ejecutar un código en C desde la terminal de Linux se utiliza el compilador `gcc` ($\textit{GNU Compiler Collection}$) mediante el siguiente comando dentro de la carpeta donde esta el archivo:

```
gcc filename.c -o filename.exe -lm
```

Esto crea el ejecutable `.exe` que puede ejecutarse llamandolo:

```
./filename.exe
```

<details>
<summary>¿Qué hace cada parte del comando?</summary>

<br>

| Parte | Significado |
|---|---|
| `gcc` | El compilador. |
| `filename.c` | El archivo fuente que se quiere compilar. |
| `-o filename.exe` | Nombre del ejecutable de salida. Sin esta opción, `gcc` lo llama `a.out`. |
| `-lm` | Enlaza la biblioteca matemática (`libm`). Es necesaria siempre que se use `math.h`: `sqrt`, `sin`, `pow`, las funciones de Bessel, etc. |

Dos opciones que vale la pena añadir mientras se aprende:

```
gcc -Wall -Wextra filename.c -o filename.exe -lm
```

`-Wall -Wextra` activan los avisos del compilador. Muchos errores clásicos de C
(variables sin inicializar, comparaciones sospechosas, valores de retorno que
faltan) aparecen ahí antes de que el programa se ejecute y dé un resultado
incorrecto en silencio.

</details>

</details>

---

<details>
<summary><b>📂 Contenido</b></summary>

<br>

| Archivo | Descripción |
|---|---|
| [`Codes/bessel_plot.c`](Codes/bessel_plot.c) | Calcula y grafica las funciones de Bessel $J_0(x)$ a $J_6(x)$. |

<details>
<summary>Codes/bessel_plot.c — detalles</summary>

<br>

Abre una tubería (`popen`) hacia `gnuplot` y le envía por ahí tanto los comandos
de configuración del gráfico como los datos, sin pasar por ningún archivo
intermedio. Usa `j0`, `j1` y `jn` de `math.h` para evaluar las siete funciones
sobre $x \in [0, 20]$.

```bash
gcc -Wall bessel_plot.c -o bessel_plot.exe -lm
./bessel_plot.exe
```

Requiere tener `gnuplot` instalado:

```bash
sudo pacman -S gnuplot     # Arch
```

</details>

</details>

