/**
 * Returns the string representation of the number. 
 * @param pValue {Float} The number to process.
 * @return {String} The string representation (example: 12856892.11111 => 12,856,892.11).
 */
var numberDisplay = function( pValue ){
	var new_val = "" ;
	if( ("" + pValue + "").indexOf(".") != -1 ){
		new_val = pValue.toFixed(2).replace(/(\d)(?=(\d{3})+\b)/g, '$1,');
	} else {
		new_val = pValue.toFixed().replace(/(\d)(?=(\d{3})+\b)/g, '$1,');
	}
	return new_val ;
}


var get_dispersion = function( values, counts ) { 
    var dispersion = new Array();
    
    // Unstack list
    unstacked_list = new Array();
    for( var idx = 0 ; idx < values.length ; idx++ ){
        for( var nb_add = 0 ; nb_add < counts[idx] ; nb_add++ ){
            unstacked_list.push( values[idx] );
        }
    }

    // Process metrics
    var nb_elt = unstacked_list.length ;
    dispersion['min'] = unstacked_list[0] ;
    dispersion['max'] = unstacked_list[nb_elt - 1];
    if( nb_elt % 2 == 0 ) {
        dispersion['median'] = unstacked_list[(nb_elt/2) -1] ;
    } else {
        dispersion['median'] = (unstacked_list[parseInt((nb_elt/2) -1)] + unstacked_list[parseInt(nb_elt/2)])/2 ;
    }
    // Deciles
    for( var idx = 1 ; idx <= 9 ; idx++ ){
        if( idx != 5 ) {
            dispersion[idx + '_decile'] = unstacked_list[Math.floor(idx*(nb_elt/10) + 0.5) -1] ;
        } else {
            dispersion['5_decile'] = dispersion['median'] ;
        }
    }
    // Quartiles
    dispersion['lower_quartile'] = unstacked_list[Math.floor((nb_elt/4) + 0.5) -1] ;
    dispersion['upper_quartile'] = unstacked_list[Math.floor((3*(nb_elt/4)) + 0.5) -1] ;
    
    return dispersion ;
};


function recreateChart(oldChart, elementId, option, theme, height = null) {
    const chartDom = document.getElementById(elementId);
    if (!chartDom) return null;

    if (oldChart) {
        oldChart.__ro?.disconnect?.(); // débrancher ResizeObserver
        oldChart.dispose();
    }

    // S'assurer que le conteneur a une taille visible
    if (!chartDom.style.height) chartDom.style.height = (height || 600) + "px";
    //if (!chartDom.style.width)  chartDom.style.width  = chartDom.clientWidth ? chartDom.clientWidth + "px" : "100%";

    // ⚡ pas de width/height fixés ici
    const chart = echarts.init(chartDom, theme, {renderer: 'canvas', devicePixelRatio: 3});
    chart.setOption(option);

    // Resize auto sur mutation du conteneur
    const ro = new ResizeObserver(() => !chart.isDisposed() && chart.resize({animation:false}));
    ro.observe(chartDom);
    chart.__ro = ro;
    return chart;
}


$('#themechoice').change(function() {
    var $select = $(this);
    var selectedIndex = $select.prop('selectedIndex');
    
    // Activer toutes les options
    $select.find('option').prop('disabled', false);

    // Désactiver l'option sélectionnée
    if (selectedIndex > 0) { // Ignorer l'option "Switch theme"
        $select.find('option').eq(selectedIndex).prop('disabled', true);
    }

    // Réinitialiser la sélection à "Switch theme"
    $select.prop('selectedIndex', 0);
});




function hexToRgba(hex, alpha = 1) {
	// Supprime le # si présent
	hex = hex.replace(/^#/, '');

	// Gestion du format court (#123 → #112233)
	if (hex.length === 3) {
		hex = hex.split('').map(c => c + c).join('');
	}

	const r = parseInt(hex.slice(0, 2), 16);
	const g = parseInt(hex.slice(2, 4), 16);
	const b = parseInt(hex.slice(4, 6), 16);

	return `rgba(${r}, ${g}, ${b}, ${alpha})`;
}


var table = function(pTitle, pCategories, pData, footer = undefined) {
    // Header
    var table_header = '';
    var table_header_line = "";

    for (var idx = 0; idx < pCategories.length; idx++) {
        table_header_line += "      <th data-sortable='true'>" + pCategories[idx] + "</th>\n";
    }

    table_header += "    <tr>\n" + table_header_line + "    </tr>\n";
    table_header = "  <thead>\n" + table_header + "  </thead>\n";
    var table_footer = footer ? "<tfoot>\n" + footer + "</tfoot>\n" : "";

    // Body
    var table_body = '';
    for (var data_idx = 0; data_idx < pData.length; data_idx++) {
        var table_body_row = "";
        for (var category_idx = 0; category_idx < pCategories.length; category_idx++) {
            if (typeof pData[data_idx][category_idx] === "number") {
                table_body_row += "      <td>" + numberDisplay(pData[data_idx][category_idx]) + "</td>\n";
            } else {
                table_body_row += "      <td>" + pData[data_idx][category_idx] + "</td>\n";
            }
        }
        table_body += "    <tr>\n" + table_body_row + "    </tr>\n";
    }
    table_body = "  <tbody>\n" + table_body + "  </tbody>\n";

    var table_caption = pTitle ? "  <caption>\n" + pTitle + "  </caption>\n" : "";

    return `
        <table
            class="table table-bordered table-striped"
            data-toggle="table"
            data-search="true"
            data-pagination="true"
            data-page-size="10"
            data-page-list='[5, 10, 20, 50, "All"]'
            data-show-export="true"
            data-export-types='["excel","csv"]'
            data-export-data-type="all"
        >
            ${table_header}
            ${table_body}
            ${table_caption}
            ${table_footer}
        </table>
    `;
};

var heatmapOption = function(data_type) {

    if (data_type == null) data_type = "clstr";

    var clean_type = {
        "clstr": "ASVs",
        "seq": "sequences"
    };

    var categories_ident = [1, 50, 80, 90, 95, 99, 100, 101];
    var categories_cover = [1, 50, 80, 90, 95, 99, 100, 101];

    var heatmap_data = get_alignment_heatmap_data(categories_ident, categories_cover, data_type)
        .map(function(item) {
            return [item[0], item[1], item[2] || 0];
        });

    //var frogsColor = style.getPropertyValue('--frogsColor').trim();

    return {
        title: {
            text: 'Number of ' + clean_type[data_type] + ' among their alignment results',
            left: 'center',
            textStyle: {fontWeight: 'normal'}
        },
        tooltip: {
            position: 'top',
            formatter: function(params) {
                return 'Identity: <b>' + get_displayed_categories(categories_ident)[params.data[0]] + '</b><br>'
                    + 'Coverage: <b>' + get_displayed_categories(categories_ident)[params.data[1]] + '</b><br>'
                    + 'Nb ' + clean_type[data_type] + ': <b>' + params.data[2] + '</b>';
            }
        },
        grid: {
            height: '70%',
            width: '70%',
            top: '15%'
        },
        xAxis: {
            type: 'category',
            data: get_displayed_categories(categories_ident),
            name: 'Identity',
            nameLocation: 'middle',
            nameGap: 30,
        },
        yAxis: {
            type: 'category',
            data: get_displayed_categories(categories_cover),
            name: 'Coverage',
            nameLocation: 'middle',
            nameGap: 50,
        },
        visualMap: {
            min: 0,
            max: Math.max(...heatmap_data.map(d => d[2])),
            calculable: false,
            orient: 'vertical',
            left: 'right',
            top: 'center',
            inRange: {
                color: ['#ffffff', getCssVar('--frogsColor')]
            },
            show: true,
            text: [
                Math.max(...heatmap_data.map(d => d[2])),
                0
            ],
            textStyle: {
                color: getCssVar('--frogsColor'),
                fontSize: 12
            }
        },
        series: [{
            name: clean_type[data_type],
            type: 'heatmap',
            data: heatmap_data,
            label: {
                show: true,
                color: '#000',
                fontSize: 12,
                formatter: function(params) {
                    return params.data[2];
                },
                textBorderColor: '#ffffff',
                textBorderWidth: 2
            },
            itemStyle: {
                borderColor: getCssVar('--frogsColor'),
                borderWidth: 1
            },
            emphasis: {
                itemStyle: {
                    shadowBlur: 10,
                    shadowColor: 'rgba(0,0,0,0.5)'
                }
            }
        }],
        toolbox: {
            feature: {
                saveAsImage: {}
            },
            right: '10%',
            top: 'top'
        }
    };
};

var histogramOption = function(pTitle, pYTitle, pCategories, pSeries, unity) {
    //var frogsColor = style.getPropertyValue('--frogsColor').trim();
    const frogsColor = getCssVar('--frogsColor');
    const frogsColor2 = getCssVar('--frogsColor2');
    return {
        title: {
            text: pTitle,
            left: 'center',
            textStyle: {fontWeight: 'normal'},
          },
          tooltip: {
            trigger: 'axis',
            axisPointer: {
              type: 'shadow'
            },
            formatter: function(params) {
              let header = `<span style="font-size:12px"><b>${params[0].axisValue}</b></span><br>`;
              let body = params.map(p => 
                `<span style="color:${p.color};">${p.seriesName}:</span> 
                 <b>${p.value} ${unity}</b><br>`
              ).join('');
              return header + body;
            }
          },
          legend: {
            top: 'bottom'
          },
          grid: {
            left: '8%',
            right: '5%',
            bottom: '10%',
            containLabel: true
          },
          //color: [frogsColor, frogsColor2],
          xAxis: {
            type: 'category',
            data: pCategories,
            axisLabel: {
              //color: frogsColor2,
              rotate: 45,
            },
            axisLine: {
              //lineStyle: { color: frogsColor2 }
            }
          },
          yAxis: {
            type: 'value',
            name: pYTitle,
            nameLocation: 'middle',
            nameGap: 40,
            axisLine: {
              //lineStyle: { color: frogsColor2 }
            },
            splitLine: {
              show: true,
              lineStyle: { color: 'rgba(0,0,0,0.1)' }
            },
            axisLabel: {
              //color: frogsColor2
            }
          },
          series: pSeries.map(serie => ({
            name: serie.name,
            type: 'bar',
            data: serie.data,
            barMaxWidth: '50%',
            emphasis: {
              focus: 'series'
            }
          })),
          toolbox: {
            feature: {
              saveAsImage: {
                title: 'Download',
                name: pTitle.replace(/\s+/g, '_')
              },
              dataZoom: {}
            },
            right: '5%',
            top: 'top'
          }
    }
}

var lineOption = function(pTitle, pXTitle, pYTitle, pXCategories, pData) {
    let xMin = Math.min(
        ...pData.flatMap(serie => serie.data.map(point => point[0]))
    );
    /*let colors = Array.from({ length: pData.length }, (_, i) =>
        `hsl(${(i * 360 / pData.length)}, 70%, 50%)`
    );*/
    return {
        title: {
            text: pTitle,
            textStyle: {fontWeight: 'normal'},
        },
        //color: ['#d87c7c', '#919e8b', '#d7ab82', '#6e7074', '#61a0a8', '#efa18d', '#787464', '#cc7e63', '#724e58', '#4b565b'],
        tooltip: {
            trigger: 'item',
            axisPointer: { show: false },
            formatter: function (params) {
                let tooltip_head = '<b>Length ' + params.value[0] + ' nt</b>';
                let tooltip_body = '<tr>' +
                '<td style="color:' + params.color + '">' + params.seriesName + ': </td>' +
                '<td>' + numberDisplay(params.value[1]) + '</td>' +
                '<td> seq</td>' +
                '</tr>';
            return tooltip_head + '<table>' + tooltip_body + '</table>';
            }
        },
        toolbox: {
            feature: {
                dataZoom: { title: { zoom: 'Zoom', back: 'Reset' } },
                saveAsImage: { title: 'Save as PNG' }
            }
        },
        xAxis: {
            type: 'value', // car on a des valeurs numériques (longueuheatmapChart_optionsrs)
            name: pXTitle,
            splitLine: {
                show: false
            },
            min:xMin,
            nameLocation: 'middle',
            //nameGap: 50,
            minInterval: 1,
            axisLabel: {
                formatter: function (value) {
                    return Math.round(value); // arrondi à l'entier le plus proche
                }
            }
        },
        yAxis: {
            type: 'value',
            name: pYTitle,
            nameLocation: 'middle',
            nameGap: 50,
            minInterval: 1,
            splitLine: {
                show: true
            },
            axisLabel: {
                formatter: function (value) {
                    return Math.round(value);
                }
            }
        },
        legend: {
            //type: 'scroll',
            type: 'plain',
            orient: 'horizontal',
            //bottom: 20,
            //height: 100,
            //pageButtonGap: 5 // espace entre les boutons de navigation
        },
        dataZoom: [
            {
                type: 'inside',   // zoom à la molette ou pinch
                xAxisIndex: 0,
                filterMode: 'filter'
            }
        ],
        series: pData.map(function (serie) {
            return {
                name: serie.name,
                type: 'line',
                data: serie.data,
                symbol: 'circle',
                symbolSize: 4,
                smooth: false,
            };
        })
    };
};

var lineOptionDualY = function(pTitle, pXTitle, x_values, y_axis_infos, my_series) {
    return {
        tooltip: {
            trigger: 'axis',
            axisPointer: { type: 'cross' },
            backgroundColor: 'rgba(255, 255, 255, 0.95)',
            borderWidth: 1,
            borderColor: '#ccc',
            textStyle: { color: '#333' },
            confine: true,
            extraCssText: 'box-shadow: 0 0 8px rgba(0,0,0,0.2); padding: 8px;',
            formatter: function (params) {
                if (!params || params.length === 0) return '';

                // Récupérer les max de chaque série
                const seqSeries = my_series.find(s => s.name === "Sequences");
                const asvSeries = my_series.find(s => s.name === "ASVs");
                const maxSeq = seqSeries.data[seqSeries.data.length - 1];
                const maxASV = asvSeries.data[asvSeries.data.length - 1];

                let tooltip = '<table style="border-collapse:collapse;">';

                params.forEach(p => {
                    const val = p.value; // juste le Y
                    let pct = 0;
                    if (p.seriesName === "Sequences") pct = maxSeq ? (val / maxSeq) * 100 : 0;
                    if (p.seriesName === "ASVs") pct = maxASV ? (val / maxASV) * 100 : 0;

                    tooltip += `
                        <tr>
                            <td style="color:${p.color};padding-right:8px;">${p.seriesName} :</td>
                            <td style="text-align:right;">
                                ${val.toLocaleString('en-US')} 
                                (${pct.toFixed(1)}%)
                            </td>
                        </tr>`;
                });

                tooltip += '</table>';
                return tooltip;
            },
            useHTML: true
        },
        title: {
            text: pTitle,
            textStyle: {fontWeight: 'normal'},
        },
        //grid: { right: '20%' },
        toolbox: {
            feature: {
                //dataView: { show: true, readOnly: false },
                dataZoom: { title: { zoom: 'Zoom', back: 'Reset' } },
                saveAsImage: { show: true }
            }
        },
        legend: {
            data: my_series.map(s => s.name)
        },
        xAxis: [
            {
                type: 'category',
                axisTick: { alignWithLabel: true },
                data: x_values
            }
        ],
        yAxis: y_axis_infos,
        series: my_series
    };
};

function boxplotOption(pTitle, pXTitle, pYTitle, pXCategories, boxplot_series) {
    return {
        title: {
            text: pTitle,
            left: 'center',
            subtext: 'N.B.: Use slider to zoom in.',
            textStyle: {
                fontWeight: 'normal'
            }
        },
        tooltip: {
            trigger: 'item',
            formatter: function (param) {
                let d = param.data;
                return [
                    `${pXCategories[param.dataIndex]}`,
                    `Min: ${d[0]}`,
                    `Q1: ${d[1]}`,
                    `Median: ${d[2]}`,
                    `Q3: ${d[3]}`,
                    `Max: ${d[4]}`
                ].join('<br/>');
            }
        },
        toolbox: {
            feature: {
                restore: {},
                saveAsImage: { title: 'Save as PNG' }
            }
        },
        xAxis: {
            type: 'category',
            name: pXTitle,
            data: pXCategories,
            boundaryGap: true,
            nameLocation: 'middle',
            nameGap: 30,
            axisPointer: {
                label: {
                    show: true,
                    backgroundColor: 'red'
                }
            }
        },
        yAxis: {
            type: 'value',
            name: pYTitle,
            min: 0,
            nameLocation: 'middle',
            nameGap: 45
        },
        grid: {
            containLabel: true,
            bottom: 0  // ajuste selon ton cas pour éviter les débordements
        },
        dataZoom: [
            /*{
                type: 'slider',
                fillerColor: "rgba(230, 234, 240, 0.4)",
                filterMode: 'none',
                yAxisIndex: 0,
                start: 0,
                end: 100,
                zoomLock: false,
                minValueSpan: 1,
                maxValueSpan: null
            },*/
            {
                type: 'slider', 
                yAxisIndex: 0, 
                zoomLock: false,
                minValueSpan: 1,
                maxValueSpan: null,
                width: 20,
                filterMode: 'none', 
                start: 0, 
                end: 100, 
                backgroundColor: "rgba(211,211,211,0.2)",
                fillerColor: "rgba(211,211,211,0.2)", 
                dataBackground: {
                      lineStyle: { color: "rgba(211,211,211,8)"},
                    areaStyle: {
                        color: "rgba(211,211,211,0.5)",
                        shadowColor: "rgba(211,211,211,0.5)"
                    }
                },
                borderColor: "rgb(211,211,211)",
                handleStyle: {
                    color: "rgba(211,211,211,0.2)"
                },
                moveHandleStyle: {
                    color: "rgba(211,211,211,1)",
                      opacity: 1
                },
                selectedDataBackground: {
                    areaStyle: {
                        color: "rgba(211,211,211,0.8)"
                    }
                },
                moveHandleSize: 4,
                emphasis: {
                    moveHandleStyle: {
                        color: "rgba(211,211,211,0.8)"
                    }
                }
            }
        ],
        series: boxplot_series.map(s => ({
            name: s.name,
            type: 'boxplot',
            boxWidth: "70%",
            data: s.data,
            /*itemStyle: {
                color: frogsColor,
                borderColor: frogsColor,
            },*/
            emphasis: {
                itemStyle: {
                    borderWidth: 2,
                    shadowBlur: 8,
                    shadowColor: 'rgba(0,0,0,0.4)'
                }
            }
        }))
    };
}

function barOption(pTitle, nb, yTitle, categories, series, unity, is_stacked) {
    const frogsColor = getCssVar('--frogsColor');
    const frogsColor2 = getCssVar('--frogsColor2');
    return {
        title: {
            text: pTitle,
            textStyle: {fontWeight: 'normal'}
        },
        tooltip: {
            trigger: 'axis',
            axisPointer: { type: 'shadow' },
            formatter: function (params) {
                let s = '<b>' + params[0].axisValue + '</b>';
                let sum = 0;
                params.forEach(function (point) {
                    s += '<br/><span style="color:' + point.color + ';">' + point.seriesName + ' : </span>'
                    + numberDisplay(point.value) + ' ' + unity;
                    if (!is_stacked) {
                        s += ' (' + (Math.round(point.value * 100 / nb * 100) / 100) + '%)';
                    }
                    sum += point.value;
                });
                if (is_stacked) {
                    s += '<br/>total : ' + numberDisplay(sum) + ' (' + (Math.round(sum * 100 / nb * 100) / 100) + '%)';
                }
                return s;
            }
        },
        legend: { show: true },
        xAxis: {
            type: 'category',
            data: categories,
            axisTick: { alignWithLabel: true }
        },
        yAxis: {
            type: 'value',
            nameLocation: 'center',
            min: 0,
            max: nb + 10,
            name: yTitle,
            splitLine: { show: true },
            axisLabel: { formatter: '{value}' }
        },
        series: series.map(s => ({
            name: s.name,
            type: 'bar',
            stack: is_stacked ? 'total' : null,
            data: s.data,
            label: {
                show: true,
                position: 'right',
                color: 'inherit',
                formatter: function (params) {
                    return numberDisplay(params.value);
                },
                fontWeight: 'bold'
            },
            // 👉 Ici on insère ton markLine
            markLine: nb ? {
                symbol: "none",
                silent: true,
                data: [{ yAxis: nb }],
                label: {
                    show: true,
                    position: "insideStartBottom",
                    padding: [0, 20, -30, -100],
                    rotate: 90,
                    color: frogsColor,
                    fontFamily: "Arial",
                    formatter: () =>
                        `Input sequences:\n${nb.toLocaleString("en-US")}`,
                },
                lineStyle: {
                    color: frogsColor,
                    type: "solid",
                    width: 1.5,
                }
            } : null
        })),
        toolbox: {
            feature: {
                saveAsImage: { title: 'Save as PNG' }
            }
        }
    };
}

function pieOption(value_1, value_2, label_1, label_2, title, unit, value_3 = null, label_3 = null) {
    const data = [
        { value: value_1, name: label_1 },
        { value: value_2, name: label_2 }
        ];
    if (value_3 !== null && label_3 !== null) {
        data.push({ value: value_3, name: label_3 });
    }
    
    let option = {
        title: {
        text: title,
        textStyle: {fontWeight: 'normal'},
        left: 'center' // 'left', 'right', 'center', ou valeur en %/px
        },
        //color: [frogsColor, frogsColor2],
        tooltip: {
        trigger: 'item' // 'item' (pour pie), 'axis' (pour bar/line)
        },
        legend: {
            show: false
        },
        toolbox: {
        feature: {
            saveAsImage: {}
        }
        },
        series: [
        {
            label: {
                color: "#000000", // ou "black", ou en hexadécimalup
                fontSize: 13,
                fontWeight: 'bold',
                fontFamily: "Arial",
                formatter: function(params) {
                    const name = params.name;
                    const value = params.value;
                    return `${name}: ${value.toLocaleString('fr-FR')}`;
                }
            },
            tooltip: {
                formatter: function (params) {
                    return `${params.name} <br>${unit}: <strong>${params.percent}%</strong>`;
                }
            },
            type: 'pie', // 'pie' est le type pour camembert
            radius: '50%', // peut être ['40%', '70%'] pour un donut
            data: data,
            itemStyle: {
                    borderColor: '#ffffff', // couleur du trait
                    borderWidth: 2 // épaisseur du trait
                    },
            emphasis: { 
                    focus: 'self',
                    blurScope: 'series',
                        itemStyle: { // paramétrage des ombres ( épaisseur), épaisseur de la bordure et de la couleur des ombres
                    borderWidth: 0, // supprime la bordure au hover
                    shadowBlur: 10,
                    shadowOffsetX: 5,
                    shadowColor: 'rgba(0, 0, 0, 0.5)'
                }

            },
            blur: {    //(opacité des différents effets de blurs)
            itemStyle: {
                opacity: 0.5
            },
            label: {
                opacity: 0.7
            }
            },
        }
        ]
    };
    return option;
}

function areaplotOption(pTitle, pXTitle, pYTitle, pXCategories, pData) {
    // Trouver le max des X
    let x_max = 0;
    for (const serie of pData) {
        for (const [x] of serie.data) {
            if (x > x_max) x_max = x;
        }
    }
    const tickInterval = Math.max(1, Math.floor(x_max / 10));

    // Créer une map des séries -> data X pour retrouver les index
    const seriesIndexMap = {};
    for (const s of pData) {
        seriesIndexMap[s.name] = s.data.map(d => d[0]);
    }

    let option = {
        title: {
            text: pTitle,
            left: 'center',
            textStyle: { fontWeight: 'normal' },
            subtext: 'N.B.: Use sliders to zoom in.'
        },
        grid: {
            left: 60,
            right: 60,
            top: 60,
            bottom: 120  // espace supplémentaire pour titre + dataZoom
        },
        tooltip: {
            trigger: 'axis',
            axisPointer: { type: 'cross' },
            useHTML: true,
            formatter: function (params) {
                if (!params?.length) return '';

                const xValue = params[0].value[0];
                let tooltip_head = `<caption><b>Clusters with size ≤ ${xValue}</b></caption>`;
                let tooltip_body = `
                    <thead><tr><th>Sequences</th><th>Clusters</th></tr></thead><tbody>
                `;

                params.forEach(p => {
                    const allX = seriesIndexMap[p.seriesName] || [];
                    const pointIndex = allX.findIndex(x => x === xValue);
                    const percCluster = (pointIndex >= 0 && allX.length > 0)
                        ? (((pointIndex + 1) / allX.length) * 100).toFixed(2)
                        : 'NA';

                    tooltip_body += `
                        <tr>
                            <td>${p.value[1].toFixed(2)}%</td>
                            <td>${percCluster}%</td>
                        </tr>
                    `;
                });

                tooltip_body += '</tbody>';
                return `<table id="tooltip-seqdepth" class="table caption-top">${tooltip_head}${tooltip_body}</table>`;
            }
        },
        toolbox: {
            feature: {
                restore: {},
                saveAsImage: { title: 'Save as PNG' },
            }
        },
        xAxis: {
            type: 'value',
            nameGap: 50,
            boundaryGap: true,
            name: pXTitle,
            nameLocation: 'middle',
            min: 1,
            max: x_max,
            interval: tickInterval
        },
        yAxis: {
            type: 'value',
            name: pYTitle,
            min: 0,
            max: 100
        },
        dataZoom: [
            {
                type: 'slider', 
                xAxisIndex: 0, 
                height: 20,
                filterMode: 'none', 
                start: 0, 
                end: 100, 
                backgroundColor: "rgba(211,211,211,0.2)",
                fillerColor: "rgba(211,211,211,0.2)", 
                dataBackground: {
                      lineStyle: { color: "rgba(211,211,211,8)"},
                    areaStyle: {
                        color: "rgba(211,211,211,0.5)",
                        shadowColor: "rgba(211,211,211,0.5)"
                    }
                },
                borderColor: "rgb(211,211,211)",
                handleStyle: {
                    color: "rgba(211,211,211,0.2)"
                },
                moveHandleStyle: {
                    color: "rgba(211,211,211,1)",
                      opacity: 1
                },
                selectedDataBackground: {
                    areaStyle: {
                        color: "rgba(211,211,211,0.8)"
                    }
                },
                moveHandleSize: 4,
                emphasis: {
                    moveHandleStyle: {
                        color: "rgba(211,211,211,0.8)"
                    }
                }
            },
            {
                type: 'slider', 
                yAxisIndex: 0, 
                width: 20,
                filterMode: 'none', 
                start: 0, 
                end: 100, 
                backgroundColor: "rgba(211,211,211,0.2)",
                fillerColor: "rgba(211,211,211,0.2)", 
                dataBackground: {
                      lineStyle: { color: "rgba(211,211,211,8)"},
                    areaStyle: {
                        color: "rgba(211,211,211,0.5)",
                        shadowColor: "rgba(211,211,211,0.5)"
                    }
                },
                borderColor: "rgb(211,211,211)",
                handleStyle: {
                    color: "rgba(211,211,211,0.2)"
                },
                moveHandleStyle: {
                    color: "rgba(211,211,211,1)",
                      opacity: 1
                },
                selectedDataBackground: {
                    areaStyle: {
                        color: "rgba(211,211,211,0.8)"
                    }
                },
                moveHandleSize: 4,
                emphasis: {
                    moveHandleStyle: {
                        color: "rgba(211,211,211,0.8)"
                    }
                }
            }
        ],
        
        series: pData.map(s => ({
            name: s.name,
            type: 'line',
            data: s.data,
            areaStyle: {},
            symbol: 'circle',
            symbolSize: 8,
            emphasis: {
                focus: 'series',
                scale: true,
                itemStyle: {
                    //borderColor: getCssVar('--frogsButtonColor'),
                    shadowColor: 'rgba(0,0,0,0.3)',
                    shadowBlur: 10,
                    shadowOffsetX: 0,
                    shadowOffsetY: 0,
                    borderWidth: 2,
                    borderColor: "#fff",
                    //color: getCssVar('--frogsColorHover')
                }
            }
        }))
    };

    return option;
}