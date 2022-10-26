# zippy.nvim

nvim plugin written in lua to help with debugging. this is inspired by the popular vscode extension [turbo console log](https://marketplace.visualstudio.com/items?itemName=ChakrounAnas.turbo-console-log).

the main functionality is to insert meaningful print statements automatically.

![image](https://user-images.githubusercontent.com/80820813/198012225-9cafe569-fa80-461e-9a58-dbbb25653f50.png)

![image](https://user-images.githubusercontent.com/80820813/198012504-11da9e43-4bdd-4982-b37f-7c8670e91d37.png)

the print statement is formatted as 
```
- file name -> line number -> variable
```

if you have [nvim navic](https://github.com/SmiteshP/nvim-navic) installed, breadcrumbs to your function are additionally added:
```
- file name -> line number -> breadcrumbs -> variable
```

## Installation

<a href="https://github.com/nvim-treesitter/nvim-treesitter" target="_blank" rel="noopener noreferrer"> Tree-sitter is required</a> <br/>

Packer:

```
use 'PatschD/zippy.nvim'
```

vim.keymap.set("n", "<leader>lg", "<cmd>lua require('zippy').insert_print()<CR>")


