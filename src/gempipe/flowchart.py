import uuid
from IPython.display import display, HTML


class Flowchart:
    def __init__(self, diagram):
        self.diagram = self.process_diagram(diagram)
        self.uid = uuid.uuid4()

    def process_diagram(self, diagram):
        diagram = diagram.replace("\n", "\\n")
        diagram = diagram.lstrip("\\n")
        diagram = diagram.replace("'", '"')
        return diagram

    def render(self, height=500, panzoom=True):
        
        panzoom_directive = 'data-zoom-on-wheel data-pan-on-drag' if panzoom else ''
        html = f"""
        <style> #outcellbox {{display: flex; justify-content: center; overflow: hidden; width: 99%; height: {height}px; background-color: #ffffff; border: 1px solid grey;}} </style>
        <script src="https://cdn.jsdelivr.net/npm/svg-pan-zoom-container@0.6.1"></script>
        <div {panzoom_directive} class="mermaid-{self.uid}" id="outcellbox"></div>  
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