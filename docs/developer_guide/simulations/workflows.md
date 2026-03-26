# Hermes workflow toolkit

The `hermesWorkflowToolkit` is the base class for simulation toolkits that support workflow pipelines:

| Feature | Description |
|---------|-------------|
| Workflow groups | Organize simulations into named groups |
| Chained steps | Chain multiple simulation steps (preprocessing → solving → post-processing) |
| Template support | Use templates for reproducible case setup |

`OFToolkit` inherits from `hermesWorkflowToolkit`, gaining workflow support automatically.
