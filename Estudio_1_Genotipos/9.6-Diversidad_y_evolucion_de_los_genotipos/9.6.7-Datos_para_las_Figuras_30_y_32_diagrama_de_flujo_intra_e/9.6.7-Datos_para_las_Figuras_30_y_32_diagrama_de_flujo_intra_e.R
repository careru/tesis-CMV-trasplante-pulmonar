# Estudio 1 · Genotipos de las glucoproteínas gB y gH de CMV
# Apartado: 9.6.7. Datos para las Figuras 30 y 32 (diagrama de flujo intra e interepisodio)
# Extracto del script completo (ver README.md del estudio)

## ---- Datos para las Figuras 30 y 32 (diagrama de flujo intra e interepisodio) ----
# Las Figuras 30 y 32 de la tesis son infografías; este diagrama recoge sus recuentos.
grViz("
digraph flujo {

  graph [layout = dot, rankdir = TB,
         fontname = 'Arial', splines = ortho,
         nodesep = 0.6, ranksep = 0.8]

  node [shape = rectangle, style = filled,
        fontname = 'Arial', fontsize = 19,
        width = 2.0, height = 0.7, margin = '0.18,0.12']

  edge [color = '#999999', arrowsize = 0.5,
        arrowhead = vee]

subgraph cluster_A {
  label = < <B>A. Variabilidad intraepisodio</B> >
  fontsize = 25
  labeljust = 'c'
  margin = 40
  color = white

    A1 [label = '234 episodios\ngenotipado gB/gH',     fillcolor = '#607D8B', fontcolor = white]
    A2 [label = '131 (56,0%)\n1 sola muestra',         fillcolor = '#90A4AE', fontcolor = white]
    A3 [label = '103 (44,0%)\n≥2 muestras',            fillcolor = '#607D8B', fontcolor = white]
    A4 [label = '61 (59,2%)\ngenotipo único',           fillcolor = '#1565C0', fontcolor = white]
    A5 [label = '42 (40,8%)\ngenotipo mixto',           fillcolor = '#BF360C', fontcolor = white]
    A6 [label = '61 (100%)\npatrón estable',            fillcolor = '#2E7D32', fontcolor = white]
    A7 [label = '7 (16,7%)\nestable',                   fillcolor = '#2E7D32', fontcolor = white]
    A8 [label = '35 (83,3%)\nfluctuante',               fillcolor = '#E65100', fontcolor = white]
    A9 [label = '34 (97,1%)\n≥1 dominante',             fillcolor = '#2E7D32', fontcolor = white]
    A10[label = '1 (2,9%)\ncambio completo',            fillcolor = '#7B1FA2', fontcolor = white]

    A1 -> A2
    A1 -> A3
    A3 -> A4
    A3 -> A5
    A4 -> A6
    A5 -> A7
    A5 -> A8
    A8 -> A9
    A8 -> A10

    {rank = same; A2; A3}
    {rank = same; A4; A5}
    {rank = same; A6; A7; A8}
    {rank = same; A9; A10}
  }

  subgraph cluster_B {
    label = < <B>B. Variabilidad interepisodio</B> >
    fontsize = 25
    labeljust = 'c'
      margin = 40
    color = white

    B1 [label = '146 pacientes\ngenotipado gB/gH',     fillcolor = '#607D8B', fontcolor = white]
    B2 [label = '90 (61,6%)\n1 solo episodio',         fillcolor = '#90A4AE', fontcolor = white]
    B3 [label = '56 (38,4%)\n>1 episodio',             fillcolor = '#607D8B', fontcolor = white]
    B4 [label = '35 (62,5%)\nmantienen patrón',        fillcolor = '#1565C0', fontcolor = white]
    B5 [label = '21 (37,5%)\ncambian patrón',          fillcolor = '#BF360C', fontcolor = white]
    B6 [label = '26 (74,3%)\nsiempre único',           fillcolor = '#1565C0', fontcolor = white]
    B7 [label = '9 (25,7%)\nsiempre mixto',            fillcolor = '#BF360C', fontcolor = white]
    B8 [label = '20 (95,2%)\n≥1 dominante',            fillcolor = '#2E7D32', fontcolor = white]
    B9 [label = '1 (4,8%)\ncambio completo',           fillcolor = '#7B1FA2', fontcolor = white]

    B1 -> B2
    B1 -> B3
    B3 -> B4
    B3 -> B5
    B4 -> B6
    B4 -> B7
    B5 -> B8
    B5 -> B9

    {rank = same; B2; B3}
    {rank = same; B4; B5}
    {rank = same; B6; B7; B8; B9}
  }
}
")


