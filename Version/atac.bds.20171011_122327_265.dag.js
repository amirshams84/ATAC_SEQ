$(function(){ // on dom ready

var cy = cytoscape({
  container: document.getElementById('cy'),
  
  style: [
    {
      selector: 'node',
      css: {
        'content': 'data(id)',
        'text-valign': 'center',
        'text-halign': 'center'
      }
    },
    {
      selector: 'edge',
      css: {
        'target-arrow-shape': 'triangle'
      }
    },
    {
      selector: ':selected',
      css: {
        'background-color': 'red',
        'line-color': 'red',
        'target-arrow-color': 'red',
        'source-arrow-color': 'red'
      }
    }
  ],
  
  elements: {
    nodes: [
      { data: { id: 'thread_Root' } },
      { data: { id: 'detect_adapter rep1', parent: 'thread_Root' } },
      { data: { id: 'detect_adapter rep1', parent: 'thread_Root' } },
      { data: { id: 'trim_adapters_PE rep1', parent: 'thread_Root' } },
      { data: { id: 'read_length rep1', parent: 'thread_Root' } },
      { data: { id: 'bowtie2_PE rep1', parent: 'thread_Root' } },
      { data: { id: 'flagstat_bam rep1', parent: 'thread_Root' } },
      { data: { id: 'dedup_bam_PE_1 rep1', parent: 'thread_Root' } },
      { data: { id: 'markdup_bam_picard rep1', parent: 'thread_Root' } },
      { data: { id: 'dedup_bam_PE_2 rep1', parent: 'thread_Root' } },
      { data: { id: 'nmsrt_bam rep1', parent: 'thread_Root' } },
      { data: { id: 'bam_to_bedpe rep1', parent: 'thread_Root' } },
      { data: { id: 'bedpe_to_tag rep1', parent: 'thread_Root' } },
      { data: { id: 'spr_PE rep1', parent: 'thread_Root' } },
      { data: { id: 'shift_tag rep1', parent: 'thread_Root' } },
      { data: { id: 'shift_tag rep1', parent: 'thread_Root' } },
      { data: { id: 'shift_tag rep1', parent: 'thread_Root' } },
      { data: { id: 'subsample_bedpe2tag rep1', parent: 'thread_Root' } },
      { data: { id: 'xcor rep1', parent: 'thread_Root' } },
      { data: { id: 'detect_adapter rep2', parent: 'thread_Root' } },
      { data: { id: 'detect_adapter rep2', parent: 'thread_Root' } },
      { data: { id: 'trim_adapters_PE rep2', parent: 'thread_Root' } },
      { data: { id: 'read_length rep2', parent: 'thread_Root' } },
      { data: { id: 'bowtie2_PE rep2', parent: 'thread_Root' } },
      { data: { id: 'flagstat_bam rep2', parent: 'thread_Root' } },
      { data: { id: 'dedup_bam_PE_1 rep2', parent: 'thread_Root' } },
      { data: { id: 'markdup_bam_picard rep2', parent: 'thread_Root' } },
      { data: { id: 'dedup_bam_PE_2 rep2', parent: 'thread_Root' } },
      { data: { id: 'nmsrt_bam rep2', parent: 'thread_Root' } },
      { data: { id: 'bam_to_bedpe rep2', parent: 'thread_Root' } },
      { data: { id: 'bedpe_to_tag rep2', parent: 'thread_Root' } },
      { data: { id: 'spr_PE rep2', parent: 'thread_Root' } },
      { data: { id: 'shift_tag rep2', parent: 'thread_Root' } },
      { data: { id: 'shift_tag rep2', parent: 'thread_Root' } },
      { data: { id: 'shift_tag rep2', parent: 'thread_Root' } },
      { data: { id: 'subsample_bedpe2tag rep2', parent: 'thread_Root' } },
      { data: { id: 'xcor rep2', parent: 'thread_Root' } },
      { data: { id: 'pool_tag pooled_rep', parent: 'thread_Root' } },
      { data: { id: 'pool_tag ppr1', parent: 'thread_Root' } },
      { data: { id: 'pool_tag ppr2', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s rep1', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s rep1-pr1', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s rep1-pr2', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s rep2', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s rep2-pr1', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s rep2-pr2', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s pooled_rep', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s ppr1', parent: 'thread_Root' } },
      { data: { id: 'macs2 n/s ppr2', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks rep1', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks rep1-pr1', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks rep1-pr2', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks rep2', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks rep2-pr1', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks rep2-pr2', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks pooled_rep', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks ppr1', parent: 'thread_Root' } },
      { data: { id: 'filt_top_peaks ppr2', parent: 'thread_Root' } },
      { data: { id: 'naive_overlap ', parent: 'thread_Root' } },
      { data: { id: 'idr2 rep1-rep2', parent: 'thread_Root' } },
      { data: { id: 'idr2 rep1-pr', parent: 'thread_Root' } },
      { data: { id: 'idr2 rep2-pr', parent: 'thread_Root' } },
      { data: { id: 'idr2 ppr', parent: 'thread_Root' } },
      { data: { id: 'FRiP rep1-rep2', parent: 'thread_Root' } },
      { data: { id: 'FRiP rep1-pr', parent: 'thread_Root' } },
      { data: { id: 'FRiP rep2-pr', parent: 'thread_Root' } },
      { data: { id: 'FRiP ppr', parent: 'thread_Root' } },
      { data: { id: 'copy file', parent: 'thread_Root' } },
      { data: { id: 'copy file', parent: 'thread_Root' } },
      { data: { id: 'copy file', parent: 'thread_Root' } },
      { data: { id: 'copy file', parent: 'thread_Root' } },
      { data: { id: 'idr final qc', parent: 'thread_Root' } },
      { data: { id: 'ataqc rep1', parent: 'thread_Root' } },
      { data: { id: 'ataqc rep2', parent: 'thread_Root' } },
      { data: { id: 'blacklist_filter peak_pooled', parent: 'thread_Root' } },
      { data: { id: 'peak_to_bigbed peak_pooled', parent: 'thread_Root' } },
      { data: { id: 'blacklist_filter peak 1', parent: 'thread_Root' } },
      { data: { id: 'peak_to_bigbed peak 1', parent: 'thread_Root' } },
      { data: { id: 'blacklist_filter peak 2', parent: 'thread_Root' } },
      { data: { id: 'peak_to_bigbed peak 2', parent: 'thread_Root' } },
      { data: { id: 'peak_to_bigbed idr peak opt', parent: 'thread_Root' } },
      { data: { id: 'peak_to_bigbed idr peak consv', parent: 'thread_Root' } },
      { data: { id: 'peak_to_bigbed peak overlap', parent: 'thread_Root' } },
      { data: { id: 'tar_all_logs', parent: 'thread_Root' } },
      { data: { id: 'peak2hammock', parent: 'thread_Root' } },
      { data: { id: 'peak2hammock', parent: 'thread_Root' } },
      { data: { id: 'peak2hammock', parent: 'thread_Root' } },
      { data: { id: 'peak2hammock', parent: 'thread_Root' } },
      { data: { id: 'peak2hammock', parent: 'thread_Root' } },
      { data: { id: 'peak2hammock', parent: 'thread_Root' } },
      { data: { id: 'report', parent: 'thread_Root' } },
      { data: { id: 'pdf2png', parent: 'thread_Root' } },
      { data: { id: 'pdf2png', parent: 'thread_Root' } },
    ],
    edges: [
      { data: { id: 'None-thread_Root', source: 'None', target: 'thread_Root' } },
      { data: { id: 'threadid-threadid', source: 'threadid', target: 'threadid' } },
      { data: { id: 'taskid-taskid', source: 'taskid', target: 'taskid' } },
    ]
  },
  
  layout: {
    name: 'breadthfirst',
  }
});

}); // on dom ready

