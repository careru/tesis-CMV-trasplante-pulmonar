# Aspectos virológicos de la infección por citomegalovirus en receptores de trasplante pulmonar

Código en R de la tesis doctoral de **Carmen Alhena Reyes Ruiz**.

- **Programa:** Doctorado en Microbiología, Universidad Autónoma de Madrid
- **Directora:** Dra. María Dolores Folgueira López
- **Centro:** Servicio de Microbiología, Hospital Universitario 12 de Octubre (Madrid)
- **Año:** 2026

## Resumen

La infección por citomegalovirus (CMV) continúa siendo una de las principales complicaciones infecciosas tras el trasplante pulmonar. Aunque la profilaxis y la monitorización han reducido la incidencia de enfermedad, persisten incertidumbres sobre el impacto clínico de la diversidad genotípica del virus y la utilidad de la secuenciación masiva para el estudio de la resistencia antiviral. Esta tesis aborda ambos aspectos para mejorar la caracterización virológica de la infección. Se realizaron dos estudios complementarios en receptores del Hospital Universitario 12 de Octubre. El primero analizó la diversidad genotípica de las glucoproteínas gB y gH mediante PCR cuantitativa genotipo-específica en 146 pacientes, 234 episodios de DNAemia y 483 muestras, evaluando su relación con la cinética viral, la evolución longitudinal de los genotipos y los desenlaces clínicos. El segundo estudió la resistencia antiviral mediante secuenciación masiva (NGS) de los genes UL97, UL54 y UL56 en 11 pacientes con DNAemia de alto nivel (≥10.000 UI/mL), analizando mutaciones de resistencia, polimorfismos naturales y variantes no descritas a distintos umbrales de detección. Las infecciones por genotipos mixtos fueron frecuentes y se asociaron con un peor control de la infección por CMV: aclaramiento viral más prolongado, mayor recurrencia, mayor número de episodios y mayor exposición viral acumulada. Los genotipos mixtos de gB se asociaron, además, con un mayor riesgo de neoplasias, especialmente de tumores sólidos no cutáneos, sin relación con otros desenlaces clásicos del trasplante. El patrón genotípico mostró una elevada estabilidad longitudinal. En el estudio de resistencia, el cribado proactivo mediante NGS identificó mutaciones de resistencia en un único paciente y numerosas variantes no descritas. En ese paciente, el análisis por debajo del umbral del 15% identificó una mutación de resistencia adicional, sin modificar el espectro de resistencia ya documentado. La diversidad genotípica del CMV constituye un biomarcador asociado a un peor control de la infección y podría contribuir a la estratificación del riesgo en el trasplante pulmonar. La NGS es una herramienta útil para monitorizar la resistencia antiviral, aunque el principal reto es interpretar la relevancia clínica de las variantes detectadas. La integración de la diversidad genotípica, la resistencia antiviral y los biomarcadores inmunológicos puede favorecer una monitorización más individualizada de estos pacientes.

**Palabras clave:** Citomegalovirus; Trasplante pulmonar; Diversidad genotípica; Resistencia antiviral.

<details>
<summary><b>Abstract (English)</b></summary>

Cytomegalovirus (CMV) infection remains one of the leading infectious complications after lung transplantation. Although prophylactic and monitoring strategies have reduced the incidence of CMV disease, important uncertainties remain regarding the clinical impact of CMV genotypic diversity and the role of next-generation sequencing (NGS) in antiviral resistance testing. This thesis addresses both aspects to improve the virological characterization of the infection. Two complementary studies were conducted in lung transplant recipients at Hospital Universitario 12 de Octubre. The first evaluated the genotypic diversity of the CMV glycoproteins gB and gH using genotype-specific quantitative PCR in a retrospective cohort of 146 patients, including 234 DNAemia episodes and 483 samples, and assessed its association with viral kinetics, longitudinal genotype dynamics, and major clinical outcomes. The second investigated antiviral resistance by NGS of the UL97, UL54, and UL56 genes in 11 patients with high-level CMV DNAemia (≥10,000 IU/mL), analyzing resistance mutations, natural polymorphisms, and previously undescribed variants at different variant-frequency thresholds. Mixed-genotype infections were common and were associated with poorer virological control, characterized by prolonged viral clearance, increased recurrence, more DNAemia episodes, and greater cumulative viral exposure. Notably, mixed gB genotypes were also associated with an increased risk of malignancy, particularly non-cutaneous solid tumors, whereas no association was observed with other classical transplant outcomes. Genotypic patterns showed high longitudinal stability. In the resistance study, proactive NGS screening identified antiviral resistance mutations in only one patient and numerous previously undescribed variants. In that patient, analysis below the 15% threshold identified an additional resistance mutation without changing the previously documented resistance spectrum. CMV genotypic diversity may represent a biomarker of poor virological control and could contribute to improved risk stratification in lung transplant recipients. NGS is a valuable tool for antiviral resistance monitoring, although the main challenge is the clinical interpretation of detected variants. Future integration of genotypic diversity, antiviral resistance, and immune biomarkers may enable more individualized monitoring of these patients.

**Keywords:** Cytomegalovirus; Lung transplant; Genotypic diversity; Antiviral resistance.

</details>

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
