import plotly
import tempfile
import plotly.tools as tls
from six import BytesIO


def export_figure(result, fs, settings, only_return_binary: bool = False):
    filename, matplot_fig, plotly_fig, df, return_data = None, None, None, None, None
    if len(result) == 5:
        filename, matplot_fig, plotly_fig, df, return_data = result
    elif len(result) == 4:
        filename, matplot_fig, plotly_fig, df = result
    elif len(result) == 2:
        filename, matplot_fig = result
    else:
        raise Exception(f'Function inside @exported_plot returned tuple of invalid size ({len(result)}).')
    if (settings.output_matplot_images and matplot_fig is not None) or (settings.create_pdf and only_return_binary and plotly_fig is None):
        bytes_buf = BytesIO()
        matplot_fig.savefig(bytes_buf, format='png')
        bytes_buf.seek(0)
        if only_return_binary:
            return bytes_buf, "png"
        with fs.open(f'{filename}.png', mode='wb') as outfile:
            outfile.write(bytes_buf.getbuffer())
        bytes_buf.close()
    if settings.output_csv and (df is not None):
        tf = tempfile.NamedTemporaryFile(suffix='.csv')
        df.to_csv(tf.name, index=False)
        tf.seek(0)
        bytes_buf = BytesIO(tf.read())
        tf.close()
        if only_return_binary:
            return bytes_buf, "csv"
        with fs.open(f'{filename}.csv', mode='wb') as outfile:
            outfile.write(bytes_buf.getbuffer())
        bytes_buf.close()
    if ((settings.create_pdf and only_return_binary) or settings.output_matplot_images)  and not settings.output_matplot_images:
        if plotly_fig is None:
            plotly_fig = tls.mpl_to_plotly(matplot_fig)
        tf = tempfile.NamedTemporaryFile(suffix='.png')
        plotly_fig.write_image(tf.name)
        tf.seek(0)
        bytes_buf = BytesIO(tf.read())
        tf.close()
        if only_return_binary:
            return bytes_buf, "png"
        with fs.open(f'{filename}.png', mode='wb') as outfile:
            outfile.write(bytes_buf.getbuffer())
        bytes_buf.close()
    if settings.output_plotly_json:
        if plotly_fig is None:
            plotly_fig = tls.mpl_to_plotly(matplot_fig)
        json_str = plotly.io.to_json(plotly_fig)
        if only_return_binary:
            return None, "json"
        with fs.open(f'{filename}.json', mode='w') as outfile:
            outfile.write(json_str)
    if settings.output_plotly_charts:
        if plotly_fig is None:
            plotly_fig = tls.mpl_to_plotly(matplot_fig)
        tf = tempfile.NamedTemporaryFile(suffix='.html')
        plotly.offline.plot(plotly_fig, filename=tf.name, auto_open=False)
        tf.seek(0)
        bytes_buf = BytesIO(tf.read())
        tf.close()
        if only_return_binary:
            return bytes_buf, "html"
        with fs.open(f'{filename}.html', mode='wb') as outfile:
            outfile.write(bytes_buf.getbuffer())
        bytes_buf.close()
    if only_return_binary:
        return None, "unknown"
    return 42