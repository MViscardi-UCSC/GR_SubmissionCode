import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

from pathlib import Path

import numpy as np
import pandas as pd

pd.set_option('display.width', 200)
pd.set_option('display.max_columns', None)

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

working_dir = Path.cwd()

comparisons = (
    ('N2', 'S5'),
    ('N2', 'S6'),
    ('S5', 'S6'),
    # ('S5', 'S6'),
)


def reprocess_pre_generated_data_for_MAPlot2(results_extra_df_dict: dict[str, pd.DataFrame], comparison_str: str,
                                             fishers_cutoff: float = 0.05, deseq2_cutoff: float = 0.05):
    results_extra_df = results_extra_df_dict[comparison_str].copy()
    comparison = comparison_str.split('|')
    comp_key = f"{comparison[0]} vs {comparison[1]}"
    results_extra_df[f'{comp_key} significant'] = ((results_extra_df[f'{comp_key} p_value'] <= fishers_cutoff)
                                                   &
                                                   (results_extra_df[f'{comp_key} p_value'] != -1))
    results_extra_df['significant'] = results_extra_df['padj'] <= deseq2_cutoff
    results_extra_df['Fishers + DESeq2 significant'] = results_extra_df['significant'] & results_extra_df[f'{comp_key} significant']
    results_extra_df['four_color_sig'] = (results_extra_df['significant'].astype(int)
                                          + results_extra_df[f'{comp_key} significant'].astype(int) * 2)
    results_extra_df['four_color_sig'] = results_extra_df['four_color_sig'].astype(str)
    results_extra_df['four_color_sig'] = results_extra_df['four_color_sig'].map({
        '0': 'black',  # Neither
        '1': 'blue',   # DESeq2
        '2': 'red',  # Fishers
        '3': 'purple',     # Both
    })
    return results_extra_df


def create_MAPlot_plotly2(res_extra_df: pd.DataFrame, comparison_str: str,
                          working_dir: Path, show=True, save_name=None, colorby='DESeq2'):
    print(f"Creating MAPlot for {comparison_str}")
    comparison = comparison_str.split('|')
    comp_key = f"{comparison[0]} vs {comparison[1]}"
    if colorby == 'DESeq2':
        color_col = 'significant'
        color_discrete_map = {True: 'red', False: 'black'}
        symbol_col = f'{comp_key} significant'
        symbol_map = {True: 'cross', False: 'circle'}
    elif colorby == 'Fishers':
        color_col = f'{comp_key} significant'
        color_discrete_map = {True: 'red', False: 'black'}
        symbol_col = 'significant'
        symbol_map = {True: 'cross', False: 'circle'}
    elif colorby == 'Both':
        color_col = 'four_color_sig'
        color_discrete_map = {'black': 'black',
                              'blue': 'blue',
                              'red': 'red',
                              'purple': 'purple'}
        symbol_col = None
        symbol_map = None
    else:
        raise ValueError(f"colorby must be 'DESeq2', 'Fishers', or 'Both', not {colorby}")
    fig = px.scatter(res_extra_df,
                     x=f"avg",
                     log_x=True,
                     y=f"log2FoldChange",
                     color=color_col,
                     color_discrete_map=color_discrete_map,
                     symbol=symbol_col,
                     symbol_map=symbol_map,
                     opacity=1.0,
                     labels={'avg': 'Average Expression',
                             'log2FoldChange': 'log<sub>2</sub>(Fold Change)',
                             'significant': 'DESeq2 Significant',
                             f'{comp_key} significant': 'Fishers Significant',
                             f'{comp_key} p_value': 'Fishers p-value',
                             'gene_name': 'Gene Name',
                             'gene_id': 'Gene ID',
                             'padj': 'Adjusted DESeq2 p-value',
                             'log10avg': 'log10(Average Expression)',
                             },
                     hover_name='gene_name',
                     hover_data=['gene_id',
                                 'log10avg',
                                 'log2FoldChange',
                                 'padj',
                                 f'{comp_key} significant',
                                 'significant',
                                 f'{comp_key} p_value',],
                     )
    fig.update_traces(marker=dict(line=dict(width=0.5,
                                            color='black')),
                      selector=dict(mode='markers'))
    fig.update_layout(
        yaxis_title='log<sub>2</sub>(Fold Change)',
        xaxis_title='Average Expression',
        title=f'<b>{comp_key}</b>',
        template="plotly_white",
    )
    if colorby == 'Both':
        new_legend_names = {
            'black': 'Neither',
            'blue': 'DESeq2',
            'red': 'Fishers',
            'purple': 'Both',
        }
        for trace in fig.data:
            trace.name = new_legend_names[trace.name]
    
    if save_name:
        save_name = save_name.rstrip('.html').rstrip('.svg').rstrip('.png')
        fig.write_html(working_dir / (save_name + '.html'))
    if show:
        fig.show()
    return fig

print("Loading initial results...")
df_dir = working_dir / 'preprocessedDataframes'
comparison_strings = [f'{comp[0]}|{comp[1]}' for comp in comparisons]
initial_results_df_dict = {}
for comp_str in comparison_strings:
    res_extra_df = pd.read_csv(df_dir / f'{comp_str}_MAPlotData.csv', index_col=0)
    initial_results_df_dict[comp_str] = res_extra_df

# Create a Dash app
app = dash.Dash(__name__, title='DESeq2 vs Fishers MAPlots')

comparison_options = [{'label': f'{comp[0]} vs {comp[1]}', 'value': f'{comp[0]}|{comp[1]}'} for comp in comparisons]
print(f"{comparison_options=}")
# Define the layout of the app
app.layout = html.Div(id='overall-app',
                      children=[
                          html.Div(id='options-div',
                                   children=[
                                       html.Div([
                                           html.Div([
                                               dcc.Markdown('**Comparison**'),
                                               dcc.Dropdown(
                                                   id='comparison-str-dropdown',
                                                   options=comparison_options,
                                                   value=comparison_options[0]['value'],
                                               ),
                                           ],
                                               style={'width': '50%'}),
                                           html.Div([
                                               dcc.Markdown('**Color By**'),
                                               dcc.RadioItems(
                                                   id='color-by-radio',
                                                   options=[
                                                       {'label': 'DESeq2', 'value': 'DESeq2'},
                                                       {'label': 'Fishers', 'value': 'Fishers'},
                                                       {'label': 'Both', 'value': 'Both'},
                                                   ],
                                                   value='DESeq2',
                                                   labelStyle={'display': 'block'}),
                                           ],
                                               style={'width': '50%'}),
                                       ],
                                           style={'display': 'flex'}
                                       ),
                                       html.Div([
                                           html.Div([
                                               dcc.Markdown('**Fisher\'s P-Value Significance Cutoff:**',
                                                            style={'margin-top': '10px'}, id='fishers-markdown'),
                                               dcc.Slider(
                                                   id='fishers-cutoff-slider',
                                                   min=0, max=0.5, value=0.05,
                                               ),
                                           ],
                                           ),
                                           html.Div([
                                               dcc.Markdown('**DESeq2\'s Adjusted P-Value Significance Cutoff:**',
                                                            style={'margin-top': '10px'}, id='deseq2-markdown'),
                                               dcc.Slider(
                                                   id='deseq2-cutoff-slider',
                                                   min=0, max=0.5, value=0.05,
                                               ),
                                           ],
                                           ),
                                       ],
                                       ),
                                   ]),
                          html.Div(id='ma-plot-div',
                                   children=[
                                       dcc.Graph(id='ma-plot'),
                                   ],
                                   ),
                      ])


# Define a callback function that will be triggered when the dropdown or slider value changes
@app.callback(
    [Output('ma-plot', 'figure'),
     Output('fishers-markdown', 'children'),
     Output('deseq2-markdown', 'children')],
    [Input('comparison-str-dropdown', 'value'),
     Input('color-by-radio', 'value'),
     Input('fishers-cutoff-slider', 'value'),
     Input('deseq2-cutoff-slider', 'value')]
)
def update_ma_plot(comparison_str,
                   color_by,
                   fishers_cutoff,
                   deseq2_cutoff):
    print(f"Updating MA Plot for {comparison_str} with fishers cutoff of {fishers_cutoff}")
    res_extra_df = reprocess_pre_generated_data_for_MAPlot2(initial_results_df_dict, comparison_str, fishers_cutoff, deseq2_cutoff)
    fig = create_MAPlot_plotly2(res_extra_df, comparison_str,
                                working_dir, show=False, colorby=color_by)
    return (fig,
            f'**Fisher\'s P-Value Significance Cutoff:** (Currently set to {fishers_cutoff})',
            f'**DESeq2\'s Adjusted P-Value Significance Cutoff:** (Currently set to {deseq2_cutoff})')


# Run the app
if __name__ == '__main__':
    # import plotly.io as pio
    # pio.renderers.default = 'browser'
    # pio.write_html(app, file='DESeq2vsFishers_MAPlots.html')
    print("Running the app!")
    app.run_server(debug=True,
                   port=8070,
                   dev_tools_ui=True,
                   dev_tools_hot_reload=True,
                   dev_tools_props_check=True,
                   )
    print("App finished running!")
