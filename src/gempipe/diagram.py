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

    def render(self):
        style = '<style> .flex-container {display: flex; justify-content: center; align-items: center;} </style>'
        html = f"""
        {style}
        <div class="mermaid-{self.uid} flex-container"></div>
        <script type="module">
            import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10.1.0/+esm'
            const graphDefinition = \'___diagram___\';
            const element = document.querySelector('.mermaid-{self.uid}');
            const {{ svg }} = await mermaid.render('graphDiv-{self.uid}', graphDefinition);
            element.innerHTML = svg;
        </script>
        """
        html = html.replace("___diagram___", self.diagram)
        return display(HTML(html))