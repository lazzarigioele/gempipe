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

    
    
    def render(self, height=500, zoom=1):
        
        html = f"""
        <style> #outcellbox {{display: flex; justify-content: center; overflow: hidden; width: 99%; height: {height}px; background-color: #ffffff; border: 1px solid grey;}} </style>
        
        <script src="https://unpkg.com/@panzoom/panzoom@4.5.1/dist/panzoom.min.js"></script>

        <div class="mermaid-{self.uid}" id="outcellbox"></div> 
        <script type="module">
            import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10.6.1/+esm';
            const graphDefinition = \'___diagram___\';
            const element = document.querySelector('.mermaid-{self.uid}');
            const {{ svg }} = await mermaid.render('graphDiv-{self.uid}', graphDefinition);
            element.innerHTML = svg;
            
            const elem = document.getElementById('graphDiv-{self.uid}');
            const panzoom = Panzoom(elem, {{maxScale: 50}});
            panzoom.zoom({zoom});
            elem.parentElement.addEventListener('wheel', panzoom.zoomWithWheel);
        </script>
        """
        html = html.replace("___diagram___", self.diagram)
        return display(HTML(html))