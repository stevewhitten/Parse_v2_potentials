class PartB {
    filename = '';
//    chartEl = null;
//    chartAucEl = null;
//    aucEl = null;
    summary = null;
//    data = null;
    downloadSummary = null;
//    downloadData = null;

    options = {
        chart: {
            showlegend: true,
            xaxis: {
                title: {
                    text: 'PS potential (Σ classifier distance P)',
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
                    title: 'Σ classifier distance P',
                    field: 'p_sum',
                    sorter: 'number',
                },
                {
                    title: 'Σ classifier distance P + U_π + U_q (trained using ∆h°)',
                    field: 'd_sum',
                    sorter: 'number',
                },
                {
                    title: 'Σ classifier distance P + U_π + U_q (trained using c_sat at 4°C)',
                    field: 'f_sum',
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
                    title: 'PS potential (Σ classifier distance P)',
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

    histogram_potential = [];

    constructor() {
//        this.chartEl = document.getElementById('chart-b');
//        this.chartAucEl = document.getElementById('chart-auc-b');
//        this.aucEl = document.getElementById('auc-b');
        this.summary = new Tabulator('#summary-b', this.options.summary);
//        this.data = new Tabulator('#data-b', this.options.data);
        this.downloadSummary = document.getElementById('download-summary-csv-b');
//        this.downloadData = document.getElementById('download-data-csv-b');

        this.downloadSummary.addEventListener('click', _ => {
            this.summary.download("csv", `${this.filename}_predicted_potentials.csv`);
        });
/*        this.downloadData.addEventListener('click', _ => {
          this.data.download("csv", `${this.filename}_percent_of_set_recall-b.csv`);
        });*/
    }

    reset() {
//        Plotly.newPlot(this.chartEl, [], this.options.chart);
//        this.aucEl.innerText = 'Recall vs human proteome accepted, AUC = 0';
        this.summary.clearData();
//        this.data.clearData();
    }

    proteinLongestRegion(protein, region) {
        var pRegions = protein.regions.filter(x => x.region == region);
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

    regionsDistSum(protein, region) {
        var residues = protein.residues.filter(x => x.region == region);
        var sum = 0;
        for (var i = 0; i < residues.length; i++) {
            sum = residues.reduce((accumulator, residue) => accumulator + residue.dist_norm, 0);
        }
        return sum;
    }

/*    processProteins(proteins) {
        this.histogram_potential.length = [];
        HU_PROT_HISTOGRAM.forEach((x, index) => {
            this.histogram_potential.push({
                PS_potential: index,
                count: 0,
                huProt_histogram_potential: x.potential,
            });
        });

        proteins.forEach(x => {
            var p_dist_sum = this.regionsDistSum(x, 'P');
            if (p_dist_sum >= HU_PROT_HISTOGRAM.length) {
                p_dist_sum = HU_PROT_HISTOGRAM.length - 1;
            }

            this.histogram_potential[Math.trunc(p_dist_sum)].count++;
        });

        for (let i = 0; i < this.histogram_potential.length; i++) {
            for (let j = i + 1; j < this.histogram_potential.length; j++) {
                this.histogram_potential[i].count += this.histogram_potential[j].count;
            }
        }
    }*/

    display(proteins, filename) {
        this.filename = filename;
     //   this.processProteins(proteins);
     //   this.displayChart(proteins);
     //   this.displayAucChart(proteins);
     //   this.displayAuc(proteins);
        this.displaySummary(proteins);
     //   this.displayData(proteins);
    }

/*    displayChart(proteins) {
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

        for (let i = 0; i < this.histogram_potential.length && i < PS_POTENTIAL_LENGTH; i++) {
            var set_ratio = 1;
            var huProt_ratio = 1;

            if (i > 0) {
                set_ratio = this.histogram_potential[i].count / proteins.length;
                huProt_ratio = this.histogram_potential[i].huProt_histogram_potential / HU_PROT_SEQUENCE_COUNT;
            }

            traces.psID.x.push(i);
            traces.psID.y.push(set_ratio * 100);

            traces.huProt.x.push(i);
            traces.huProt.y.push(huProt_ratio * 100);
        }

        const config = {
            toImageButtonOptions: {
                filename: `B`,
            },
        };

        const data = [traces.psID, traces.huProt];
        Plotly.newPlot(this.chartEl, data, this.options.chart, config);
    }*/
/*    displayAucChart(proteins) {
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

        for (let i = 0; i < this.histogram_potential.length; i++) {
            var set_ratio = 1;
            var huProt_ratio = 1;

            if (i > 0) {
                set_ratio = this.histogram_potential[i].count / proteins.length;
                huProt_ratio = this.histogram_potential[i].huProt_histogram_potential / HU_PROT_SEQUENCE_COUNT;
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
    }*/
/*    displayAuc(proteins) {
        var auc = 0;

        for (let i = 0; i < this.histogram_potential.length - 1; i++) {
            var x1 = this.histogram_potential[i].huProt_histogram_potential / HU_PROT_SEQUENCE_COUNT;
            var x2 = this.histogram_potential[i + 1].huProt_histogram_potential / HU_PROT_SEQUENCE_COUNT;
            var y1 = this.histogram_potential[i].count / proteins.length;
            var y2 = this.histogram_potential[i + 1].count / proteins.length;

            auc += (x1 - x2) * (y1 + y2) / 2.0;
        }

        this.aucEl.innerText = `Recall plot AUC = ${auc.toFixed(4)}`;
    }*/
    displaySummary(proteins) {
        var data = [];

        proteins.forEach((x, i) => {
            const longest = {
                P: {
                    length: 0,
                    start: 0,
                    end: 0,
                },
                D: {
                    length: 0,
                    start: 0,
                    end: 0,
                },
                F: {
                    length: 0,
                    start: 0,
                    end: 0,
                },
            };

            const longestPRegion = this.proteinLongestRegion(x, 'P');
            const longestPRegionLength = this.regionLength(longestPRegion);
            const longestDRegion = this.proteinLongestRegion(x, 'D');
            const longestDRegionLength = this.regionLength(longestDRegion);
            const longestFRegion = this.proteinLongestRegion(x, 'F');
            const longestFRegionLength = this.regionLength(longestFRegion);

            if (longestPRegionLength >= 20) {
                longest.P.length = longestPRegionLength;
                longest.P.start = longestPRegion.start + 1; // 0 index
                longest.P.end = longestPRegion.end;
            }
            if (longestDRegionLength >= 20) {
                longest.D.length = longestDRegionLength;
                longest.D.start = longestDRegion.start + 1; // 0 index
                longest.D.end = longestDRegion.end;
            }
            if (longestFRegionLength >= 20) {
                longest.F.length = longestFRegionLength;
                longest.F.start = longestFRegion.start + 1; // 0 index
                longest.F.end = longestFRegion.end;
            }

            data.push({
                index: (i + 1),
                length: x.sequence.length,
                p_sum: x.residues.filter(x => x.region == 'P').reduce((a, c) => a + c.dist_norm, 0).toFixed(1),
                longest_ps_idr: longest.P.length,
                p_start: longest.P.start,
                p_end: longest.P.end,
                d_sum: x.residues.filter(x => x.region_pi_q == 'P').reduce((a, c) => a + c.dist_norm_pi_q, 0).toFixed(1),
                longest_nonPS_IDR: longest.D.length,
                d_start: longest.D.start,
                d_end: longest.D.end,
                f_sum: x.residues.filter(x => x.region_pi_q_csat == 'P').reduce((a, c) => a + c.dist_norm_pi_q_csat, 0).toFixed(1),
                longest_Fdomain: longest.F.length,
                f_start: longest.F.start,
                f_end: longest.F.end,
                uniprot_id: x.idData.UnitProtId,
                gene_id: x.idData.GeneID,
                name: x.idData.Protein,
            });
        });

        this.summary.setData(data);
    }
/*    displayData(proteins) {
        const data = [];

        for (let i = 0; i < this.histogram_potential.length; i++) {
            var set_ratio = 1;
            var huProt_ratio = 1;

            if (i > 0) {
                set_ratio = this.histogram_potential[i].count / proteins.length;
                huProt_ratio = this.histogram_potential[i].huProt_histogram_potential / HU_PROT_SEQUENCE_COUNT;
            }

            data.push({
                index: i,
                count: this.histogram_potential[i].count,
                percent_of_set: set_ratio * 100,
                percent_of_set_huprot: huProt_ratio * 100,
                huprot_accepted: huProt_ratio,
                recall: set_ratio,
            });
        }

        this.data.setData(data);
    }*/
}
