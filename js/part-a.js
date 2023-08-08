class PartA {
    filename = '';
    chartEl = null;
    chartAucEl = null;
    aucEl = null;
    summary = null;
    data = null;
    downloadSummary = null;
    downloadData = null;

    options = {
        chart: {
            showlegend: true,
            xaxis: {
                title: {
                    text: 'length of longest predicted PS IDR',
                },
            },
            yaxis: {
                title: {
                    text: '% of set',
                },
            },
            shapes: [],
        },
        aucChart: {
            showlegend: true,
            xaxis: {
                title: {
                    text: 'human proteome accepted',
                },
            },
            yaxis: {
                title: {
                    text: 'recall',
                },
            },
            shapes: [],
        },
        summary: {
            autoResize: true,
            maxHeight: '300px',
            layout: 'fitColumns',
            columns: [
                {
                    title: 'Sequence #',
                    field: 'index',
                    sorter: 'number',
                },
                {
                    title: 'Sequence length',
                    field: 'length',
                    sorter: 'number',
                },
                {
                    title: 'length of longest PS IDR',
                    field: 'longest_ps_idr',
                    sorter: 'number',
                },
                {
                    title: 'first residue of longest PS IDR',
                    field: 'start',
                    sorter: 'number',
                },
                {
                    title: 'last residue of longest PS IDR',
                    field: 'end',
                    sorter: 'number',
                },
                {
                    title: 'UniProtKB accession ID',
                    field: 'uniprot_id',
                    sorter: 'string',
                },
                {
                    title: 'Gene',
                    field: 'gene_id',
                    sorter: 'string',
                },
                {
                    title: 'Protein',
                    field: 'name',
                    sorter: 'string',
                },
            ],
        },
        data: {
            autoResize: true,
            maxHeight: '300px',
            layout: 'fitColumns',
            columns: [
                {
                    title: 'length of longest PS IDR',
                    field: 'index',
                    sorter: 'number',
                },
                {
                    title: 'count',
                    field: 'count',
                    sorter: 'number',
                },
                {
                    title: '% of set, input file',
                    field: 'percent_of_set',
                    sorter: 'number',
                },
                {
                    title: '% of set, human proteome',
                    field: 'percent_of_set_huprot',
                    sorter: 'number',
                },
                {
                    title: 'human proteome accepted',
                    field: 'huprot_accepted',
                    sorter: 'number',
                },
                {
                    title: 'input file accepted (recall)',
                    field: 'recall',
                    sorter: 'string',
                },
            ],
        }
    };

    histogram_length = [];

    constructor() {
        this.chartEl = document.getElementById('chart-a');
        this.chartAucEl = document.getElementById('chart-auc-a');
        this.aucEl = document.getElementById('auc-a');
        //this.summary = new Tabulator('#summary-a', this.options.summary);
        this.data = new Tabulator('#data-a', this.options.data);
        //this.downloadSummary = document.getElementById('download-summary-csv-a');
        this.downloadData = document.getElementById('download-data-csv-a');

        //this.downloadSummary.addEventListener('click', _ => {
        //    this.summary.download("csv", `${this.filename}_summary-a.csv`);
        //});
        this.downloadData.addEventListener('click', _ => {
            this.data.download("csv", `${this.filename}_percent_of_set_recall-a.csv`);
        });
    }

    reset() {
        Plotly.newPlot(this.chartEl, [], this.options.chart);
        this.aucEl.innerText = 'Recall vs human proteome accepted, AUC = 0';
        this.summary.clearData();
        this.data.clearData();
    }

    proteinLongestRegion(protein) {
        var pRegions = protein.regions.filter(x => x.region == 'P');
        var longestRegion = null;
        if (pRegions.length > 0) {
            longestRegion = pRegions.reduce((prev, curr) => (prev.end - prev.start) > (curr.end - curr.start) ? prev : curr);
        }
        return longestRegion;
    }
    regionLength(region) {
        if (region == null) {
            return 0;
        }
        return region.end - region.start;
    }

    display(proteins, filename) {
        this.filename = filename;
        this.processProteins(proteins);
        this.displayChart(proteins);
        this.displayAucChart(proteins);
        this.displayAuc(proteins);
        this.displaySummary(proteins);
        this.displayData(proteins);
    }

    processProteins(proteins) {
        this.histogram_length.length = [];
        HU_PROT_HISTOGRAM.forEach((x, index) => {
            this.histogram_length.push({
                longest_PS_IDR: index,
                count: 0,
                huProt_histogram_length: x.length,
            });
        });

        proteins.forEach(x => {
            var longestRegionValue = this.regionLength(this.proteinLongestRegion(x));
            if (longestRegionValue >= HU_PROT_HISTOGRAM.length) {
                longestRegionValue = HU_PROT_HISTOGRAM.length - 1;
            }

            this.histogram_length[longestRegionValue].count++;
        });

        for (let i = 0; i < this.histogram_length.length; i++) {
            for (let j = i + 1; j < this.histogram_length.length; j++) {
                this.histogram_length[i].count += this.histogram_length[j].count;
            }
        }
    }

    displayChart(proteins) {
        const traces = {
            psID: {
                name: 'input file',
                type: 'scatter',
                mode: 'lines',
                hoverinfo: 'x+y',
                line: {
                    color: '#50daff',
                },
                x: [],
                y: [],
            },
            huProt: {
                name: 'human proteome',
                type: 'scatter',
                mode: 'lines',
                hoverinfo: 'x+y',
                line: {
                    color: '#262626',
                },
                x: [],
                y: [],
            },
        };

        for (let i = 0; i < this.histogram_length.length && i < PREDICTED_PS_LENGTH; i++) {
            var set_ratio = 1;
            var huProt_ratio = 1;

            if (i > 0) {
                set_ratio = this.histogram_length[i].count / proteins.length;
                huProt_ratio = this.histogram_length[i].huProt_histogram_length / HU_PROT_SEQUENCE_COUNT;
            }

            traces.psID.x.push(i);
            traces.psID.y.push(set_ratio * 100);

            traces.huProt.x.push(i);
            traces.huProt.y.push(huProt_ratio * 100);
        }

        const config = {
            toImageButtonOptions: {
                filename: `% of set to predidcted PS region length`,
            },
        };

        const data = [traces.psID, traces.huProt];
        Plotly.newPlot(this.chartEl, data, this.options.chart, config);
    }
    displayAucChart(proteins) {
        const traces = {
            psID: {
                name: 'input file',
                type: 'scatter',
                mode: 'lines',
                hoverinfo: 'x+y',
                line: {
                    color: '#50daff',
                },
                x: [],
                y: [],
            },
            huProt: {
                name: 'human proteome',
                type: 'scatter',
                mode: 'lines',
                hoverinfo: 'x+y',
                line: {
                    color: '#262626',
                },
                x: [],
                y: [],
            },
        };

        for (let i = 0; i < this.histogram_length.length; i++) {
            var set_ratio = 1;
            var huProt_ratio = 1;

            if (i > 0) {
                set_ratio = this.histogram_length[i].count / proteins.length;
                huProt_ratio = this.histogram_length[i].huProt_histogram_length / HU_PROT_SEQUENCE_COUNT;
            }

            traces.psID.x.push(huProt_ratio);
            traces.psID.y.push(set_ratio);

            traces.huProt.x.push(huProt_ratio);
            traces.huProt.y.push(huProt_ratio);
        }

        const config = {
            toImageButtonOptions: {
                filename: `B-C`,
            },
        };

        const data = [traces.psID, traces.huProt];
        Plotly.newPlot(this.chartAucEl, data, this.options.aucChart, config);
    }
    displayAuc(proteins) {
        var auc = 0;

        for (let i = 0; i < this.histogram_length.length - 1; i++) {
            var x1 = this.histogram_length[i].huProt_histogram_length / HU_PROT_SEQUENCE_COUNT;
            var x2 = this.histogram_length[i + 1].huProt_histogram_length / HU_PROT_SEQUENCE_COUNT;
            var y1 = this.histogram_length[i].count / proteins.length;
            var y2 = this.histogram_length[i + 1].count / proteins.length;

            auc += (x1 - x2) * (y1 + y2) / 2.0;
        }

        this.aucEl.innerText = `Recall plot AUC = ${auc.toFixed(4)}`;
    }
    displaySummary(proteins) {
        var data = [];

        proteins.forEach((x, i) => {
            const longest_ps_idr = {
                length: 0,
                start: 0,
                end: 0,
            };

            const longestRegion = this.proteinLongestRegion(x);
            const longestRegionLength = this.regionLength(longestRegion);

            if (longestRegionLength >= 20) {
                longest_ps_idr.length = longestRegionLength;
                longest_ps_idr.start = longestRegion.start + 1; // 0 index
                longest_ps_idr.end = longestRegion.end;
            }

            data.push({
                index: (i + 1),
                length: x.sequence.length,
                longest_ps_idr: longest_ps_idr.length,
                start: longest_ps_idr.start,
                end: longest_ps_idr.end,
                uniprot_id: x.idData.UnitProtId,
                gene_id: x.idData.GeneID,
                name: x.idData.Protein,
            });
        });

        //this.summary.setData(data);
    }
    displayData(proteins) {
        const data = [];

        for (let i = 0; i < this.histogram_length.length; i++) {
            var set_ratio = 1;
            var huProt_ratio = 1;

            if (i > 0) {
                set_ratio = this.histogram_length[i].count / proteins.length;
                huProt_ratio = this.histogram_length[i].huProt_histogram_length / HU_PROT_SEQUENCE_COUNT;
            }

            data.push({
                index: i,
                count: this.histogram_length[i].count,
                percent_of_set: set_ratio * 100,
                percent_of_set_huprot: huProt_ratio * 100,
                huprot_accepted: huProt_ratio,
                recall: set_ratio,
            });
        }

        this.data.setData(data);
    }
}
