# IonFlow for Galaxy #

Galaxy tool for processing and analysis of ionomics data , based on the
modification of R package [IonFlow](https://github.com/AlinaPeluso/MetaboFlow). 

## Installation ##

- Install [Galaxy](https://github.com/galaxyproject/galaxy) under Linux.

- Use `git` to clone this tool

  ```bash
  git clone https://github.com/wanchanglin/ionflow.git
  ```

- Add this tool's location into Galaxy' tool config file:
  `~/Galaxy/config/tool_conf.xml`. For example, one simplified
  `tool_conf.xml` looks like:

  ```xml
  <?xml version='1.0' encoding='utf-8'?>
  <toolbox monitor="true">
    
    <section id="getext" name="Get Data">
      <tool file="data_source/upload.xml" />
    </section>
    
    <section id="MyTools" name="My Tools">
      <tool file="/path/to/ionflow/ionflow.xml" />
    </section>

  </toolbox>
  ```

- Test data are in `test-data`.

## Authors, contributors & contacts ##

- Wanchang Lin (w.lin@imperial.ac.uk), Imperial College London
- Jacopo Iacovacci (j.iacovacci@imperial.ac.uk), Imperial College London
- Robert Glen (r.glen@imperial.ac.uk), Imperial College London
- Julian Griffin (julian.griffin@imperial.ac.uk), Imperial College London
