import uuid
from IPython.display import display, HTML


class Pipeflow:
    def __init__(self, diagram):
        self.diagram = self.process_diagram(diagram)
        self.uid = uuid.uuid4()

    def process_diagram(self, diagram):
        diagram = diagram.replace("\n", "\\n")
        diagram = diagram.lstrip("\\n")
        diagram = diagram.replace("'", '"')
        return diagram

    def render(self, height=500, panzoom='true'):
        html = f"""
        <style> #outcellbox {{display: flex; justify-content: center; width: 100%; height: {height}px; background-color: #ffffff;}} </style>
        <style> #outcellbox svg {{width: 100%; height: 100%;}} </style>
        <div class="mermaid-{self.uid}" id="outcellbox"></div>
        <script src='https://unpkg.com/panzoom@9.4.0/dist/panzoom.min.js' query='#graphDiv-{self.uid}' name='pz'></script>
        <script type="module">
            import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10.6.1/+esm';
            const graphDefinition = \'___diagram___\';
            const element = document.querySelector('.mermaid-{self.uid}');
            const {{ svg }} = await mermaid.render('graphDiv-{self.uid}', graphDefinition);
            element.innerHTML = svg;
        </script>
        """
        html = html.replace("___diagram___", self.diagram)
        return display(HTML(html))