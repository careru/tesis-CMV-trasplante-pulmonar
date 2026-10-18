# Aspectos virológicos de la infección por citomegalovirus en receptores de trasplante pulmonar

Código en R de la tesis doctoral de **Carmen Alhena Reyes Ruiz**.

- **Programa:** Doctorado en Microbiología, Universidad Autónoma de Madrid
- **Directora:** Dra. María Dolores Folgueira López
- **Centro:** Servicio de Microbiología, Hospital Universitario 12 de Octubre (Madrid)
- **Año:** 2026

> [!IMPORTANT]
> **Tesis en proceso de depósito, pendiente de defensa.** Este repositorio se pone a
> disposición de los revisores y del tribunal para la consulta del código de análisis.
> Los resultados no son definitivos hasta la defensa y la publicación de la tesis, y
> pueden estar sujetos a correcciones. No se autoriza su cita, reproducción ni
> reutilización sin el permiso expreso de la autora.

## Estructura del repositorio

| Estudio | Contenido |
|---|---|
| [Estudio 1 · Genotipos](Estudio_1_Genotipos/README.md) | Genotipado de gB y gH, cinética viral y desenlaces clínicos (Resultados, apartado 9) |
| [Estudio 2 · NGS](Estudio_2_NGS/README.md) | Variantes de resistencia en UL97, UL54 y UL56 mediante NGS (Resultados, apartado 10) |

Dentro de cada estudio hay una carpeta por apartado de Resultados (9.1, 9.2... / 10.1, 10.2)
y, dentro de ella, un archivo de código por figura o tabla. El script completo de cada
estudio está en la raíz de su carpeta.

## Protección de datos

**Qué contiene el repositorio.** Solo código en R. No se ha subido ningún archivo de datos:
ni las bases clínicas, ni los resultados de genotipado, ni las cargas virales, ni los ficheros
de secuenciación (NGS), ni las tablas o figuras generadas. El archivo `.gitignore` impide
además que se suban por error archivos `.xlsx`, `.xls`, `.dta`, `.csv`, `.RData` o `.rds`.

**Identificadores.** El código original contenía números de historia clínica (NHC) para
seleccionar pacientes. En esta copia se han sustituido por seudónimos (`NHC_RTP-1`,
`NHC_anonimizado`...) y se han suprimido las fechas personales (fecha de trasplante).
La correspondencia entre seudónimos e historias clínicas no está en el repositorio: se
conserva únicamente en el Hospital Universitario 12 de Octubre. Antes de cada subida, un
control automático busca cualquier NHC original en el código y, si encuentra alguno,
cancela la publicación.

**Información clínica presente en el código.** Solo valores agregados ya publicados en la
tesis (recuentos y umbrales) y, para la paciente RTP-2, la evolución de la carga viral y de
los tratamientos expresada en días postrasplante, tal como aparece en las Figuras 49 y 50.

**Marco ético y legal.** El estudio fue aprobado por el Comité de Ética de la Investigación
con Medicamentos del Hospital Universitario 12 de Octubre (CEIm 19/490), con exención de
consentimiento informado por su carácter retrospectivo. Los datos se trataron de forma
confidencial y seudonimizada, conforme al Reglamento General de Protección de Datos
(RGPD, UE 2016/679) y a la Ley Orgánica 3/2018 de Protección de Datos Personales y garantía
de los derechos digitales.

**Uso del código.** Los extractos sirven para consultar cómo se hizo cada análisis, pero no
se ejecutan sin los datos originales, que no son públicos. Para cualquier consulta sobre
los datos, contactar con la autora.
